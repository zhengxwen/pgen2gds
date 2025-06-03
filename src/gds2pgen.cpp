// ===========================================================
//
// gds2pgen.cpp: Format Conversion from PLINK2 PGEN to GDS
//
// Copyright (C) 2025    Xiuwen Zheng (zhengx@u.washington.edu)
//
// gds2pgen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// gds2pgen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with gds2pgen.
// If not, see <http://www.gnu.org/licenses/>.


#include <R_GDS_CPP.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define class xclass
#define private xprivate
#include <R_ext/Connections.h>
#undef class
#undef private

#include <vector>

using namespace std;
using namespace CoreArray;


// ===================================================================== //

/// Progress object
class COREARRAY_DLL_LOCAL CProgress
{
public:
	CProgress();
	CProgress(C_Int64 count, SEXP conn);

	void Reset(C_Int64 count);
	void Forward(C_Int64 val);
	void ShowProgress();

protected:
	C_Int64 fTotalCount;  ///< the total number
	C_Int64 fCounter;  ///< the current counter
	double _start, _step;
	C_Int64 _hit;
	vector< pair<double, time_t> > _timer;
	time_t _start_time, _last_time, _last_check_time;
	Rconnection progress_conn;
};

static const int PROGRESS_BAR_CHAR_NUM = 50;
static const double S_MIN  =  60;
static const double S_HOUR =  60 * S_MIN;
static const double S_DAY  =  24 * S_HOUR;
static const double S_YEAR = 365 * S_DAY;

static const char *time_str(double s)
{
	if (R_FINITE(s))
	{
		static char buffer[64];
		if (s < S_MIN)
			snprintf(buffer, sizeof(buffer), "%.0fs", s);
		else if (s < S_HOUR)
			snprintf(buffer, sizeof(buffer), "%.1fm", s/S_MIN);
		else if (s < S_DAY)
			snprintf(buffer, sizeof(buffer), "%.1fh", s/S_HOUR);
		else if (s < S_YEAR)
			snprintf(buffer, sizeof(buffer), "%.1fd", s/S_DAY);
		else
			snprintf(buffer, sizeof(buffer), "%.1f years", s/S_YEAR);
		return buffer;
	} else
		return "---";
}

extern "C"
{
	static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }
}

static bool CheckInterrupt()
{
	// this will call the above in a top-level context so it won't longjmp-out of your context
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

inline static void put_text(Rconnection conn, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	(*conn->vfprintf)(conn, fmt, args);
	va_end(args);
}



CProgress::CProgress()
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = NULL;
}

CProgress::CProgress(C_Int64 count, SEXP conn)
{
	fTotalCount = 0;
	fCounter = 0;
	progress_conn = (!conn || Rf_isNull(conn)) ? NULL : R_GetConnection(conn);
	Reset(count);
}

void CProgress::Reset(C_Int64 count)
{
	bool flag = (fTotalCount==0) || (fCounter > 0);
	fTotalCount = count;
	fCounter = 0;
	if (count > 0)
	{
		int n = 100;
		if (n > count) n = count;
		if (n < 1) n = 1;
		_start = _step = (double)count / n;
		_hit = (C_Int64)(_start);
		double percent = (double)fCounter / count;
		time_t s; time(&s);
		_start_time = _last_time = _last_check_time = s;
		_timer.clear();
		_timer.reserve(128);
		_timer.push_back(pair<double, time_t>(percent, s));
		if (flag) ShowProgress();
	}
}

void CProgress::Forward(C_Int64 val)
{
	if (fTotalCount > 0)
	{
		fCounter += val;
		if (fCounter >= _hit)
		{
			do {
				_start += _step;
				_hit = (C_Int64)(_start);
			} while (fCounter >= _hit);
			if (_hit > fTotalCount) _hit = fTotalCount;
			ShowProgress();
		}
		// check whether user interrupts
		time_t now; time(&now);
		if (difftime(now, _last_check_time) > 0.25)
		{
			_last_check_time = now;
			if (CheckInterrupt())
				throw ErrCoreArray("User interrupts!");
		}
	}
}

void CProgress::ShowProgress()
{
	if (fTotalCount > 0)
	{
		char bar[PROGRESS_BAR_CHAR_NUM + 1];
		double p = (double)fCounter / fTotalCount;
		int n = (int)round(p * PROGRESS_BAR_CHAR_NUM);
		memset(bar, '.', sizeof(bar));
		memset(bar, '=', n);
		if ((fCounter > 0) && (n < PROGRESS_BAR_CHAR_NUM))
			bar[n] = '>';
		bar[PROGRESS_BAR_CHAR_NUM] = 0;

		// ETC: estimated time to complete
		n = (int)_timer.size() - 20;  // 20% as a sliding window size
		if (n < 0) n = 0;

		time_t now; time(&now);
		_timer.push_back(pair<double, time_t>(p, now));

		// in seconds
		double interval = difftime(now, _last_time);
		double s = difftime(now, _timer[n].second);
		double diff = p - _timer[n].first;
		if (diff > 0)
			s = s / diff * (1 - p);
		else
			s = R_NaN;
		p *= 100;

		// show
		_last_time = now;
		if (fCounter >= fTotalCount)
		{
			s = difftime(_last_time, _start_time);
			if (!progress_conn)
				Rprintf("\r[%s] 100%%, completed in %s\n", bar, time_str(s));
			else {
				put_text(progress_conn, "[%s] 100%%, completed in %s\n", bar, time_str(s));
				(*progress_conn->fflush)(progress_conn);
			}
		} else if ((interval >= 5) || (fCounter <= 0))
		{
			if (!progress_conn)
				Rprintf("\r[%s] %2.0f%%, ETC: %s        ", bar, p, time_str(s));
			else {
				put_text(progress_conn, "[%s] %2.0f%%, ETC: %s\n", bar, p, time_str(s));
				(*progress_conn->fflush)(progress_conn);
			}
			// fflush(stdout);
		}
	}
}



extern "C"
{

/// Import a pgen file
COREARRAY_DLL_EXPORT SEXP SEQ_PGEN_Allele_Import(
	SEXP R_read_fc, SEXP R_allele_buf, SEXP R_ii, SEXP R_phase_buf, SEXP R_env,
	SEXP gds_root, SEXP Start, SEXP Count, SEXP progfile, SEXP Verbose)
{
	int start = Rf_asInteger(Start);
	if (start < 1) start = 1;
	int count = Rf_asInteger(Count);
	const bool verbose = Rf_asLogical(Verbose)==TRUE;
	const C_UInt8 ONE = 1;

	COREARRAY_TRY

		// GDS nodes
		PdAbstractArray Root = GDS_R_SEXP2Obj(gds_root, FALSE);
		PdAbstractArray varGeno = GDS_Node_Path(Root, "genotype/data", FALSE);
		PdAbstractArray varGenoLen = GDS_Node_Path(Root, "genotype/@data", FALSE);
		PdAbstractArray varPhase = GDS_Node_Path(Root, "phase/data", FALSE);

		const int *pdim = INTEGER(GET_DIM(R_allele_buf));
		const int num_ploidy = pdim[0], num_sample = pdim[1];
		const size_t ntot = num_ploidy * num_sample;
		int *ptr_allele_buf = INTEGER(R_allele_buf);
		int *ptr_phase_buf = INTEGER(R_phase_buf);

		// progress information
		CProgress prog((verbose || !Rf_isNull(progfile)) ? count : -1, progfile);

		// for-loop
		for (size_t idx=start; count > 0; idx++, count--)
		{
			// call ReadAlleles
			INTEGER(R_ii)[0] = idx;
			Rf_eval(R_read_fc, R_env);
			// replace missing value
			for (size_t i=0; i < ntot; i++)
			{
				if (ptr_allele_buf[i] == NA_INTEGER)
					ptr_allele_buf[i] = 3;
			}
			// append genotypes
			GDS_Array_AppendData(varGeno, ntot, ptr_allele_buf, svInt32);
			GDS_Array_AppendData(varGenoLen, 1, &ONE, svUInt8);
			GDS_Array_AppendData(varPhase, num_sample, ptr_phase_buf, svInt32);
			// update progress
			prog.Forward(1);
		}

	COREARRAY_CATCH
}


/// initialize the package
COREARRAY_DLL_EXPORT void R_init_gds2pgen(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }
	static R_CallMethodDef callMethods[] =
	{
		CALL(SEQ_PGEN_Allele_Import, 10),
		{ NULL, NULL, 0 }
	};
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	Init_GDS_Routines();
}

} // extern "C"
