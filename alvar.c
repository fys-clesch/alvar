/**
 * Author: Clemens Sch\"afermeier
 * Release: 2014/11/20
 * Contact: clesch(at)fysik.dtu.dk
 * Purpose: Load a data file, compute the Allan variance of each odd column and print the results to a file.
 *
 * This file is part of alvar.
 * 
 * alvar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * alvar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with alvar.  If not, see <http://www.gnu.org/licenses/>.
 */

/**@TODO:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <omp.h>

#ifdef __unix__
 #ifndef __USE_BSD /**< for S_IFDIR */
  #define __USE_BSD
 #endif
#elif defined(__WIN32__) || defined(_WIN32) || defined(__MSDOS__)
#endif

#include <sys/stat.h>

typedef unsigned int uint;
typedef unsigned char uchar;

#define ERR_ARG __LINE__, __func__
#define ERR_INTRO fprintf(stderr, "DIRE STRAITS in line %i from '%s': "

/* Internal header */

uchar parse_args(const int argc, char **argv,
				char *__restrict fname_in, char *__restrict fname_out,
				const uint len_fname, uchar *__restrict time_data);
double *get_alvar(const double *__restrict x, const uint len, uint *__restrict n_tau);
double *get_alvar_timeseries(const double *__restrict x, const uint len, uint *__restrict len_res, const double tau_0);
double get_sum_sqdev(const double *__restrict x, const uint tau, const uint len);
uint fpfile2double(const char *__restrict target, double *__restrict m, const uint row, const uint col, const uchar chatty);
void count_entries(char *__restrict fname,
				uint *__restrict nrow, uint *__restrict ncol,
				uint *__restrict cmt_lines);
uchar isfile(const char *szfile);
void double2fpfile(const char *__restrict fname, const double *__restrict out, const uint row, const uint col);
double *malloc_double(const uint size);
void show_progressbar(const uint inc, const uchar order, const uint ntot);
void show_progress(const uint n, const uchar order, const uint ntot);
char *replace_char(const char *str_in, const char find, const char replace);
void plotdata(const char *__restrict fname, const uint ncol);

/* Definitions */

int main(int argc, char **argv)
{
	const uint s_fname = 128;
	uchar use_timeseries_def = 1;
	char fname_in[s_fname], fname_out[s_fname];
	if(!parse_args(argc, argv, fname_in, fname_out, s_fname, &use_timeseries_def))
		return EXIT_FAILURE;

	if(!strncmp(fname_in, "NULL", s_fname)) /**> Debugging only */
	{
		double dummy[599];
		uint len = 599, tau, i;

		for(i = 0; i < len; i++)
		{
			dummy[i] = rand();
			fprintf(stdout, "input[%u]: %g\n", i, dummy[i]);
		}

		double *alvar = get_alvar(dummy, len, &tau);

		for(i = 0; i < tau; i++)
			fprintf(stdout, "alvar[%u]: %g\n", i, log(alvar[i]));

		free(alvar);
	}
	else /**> Normal operation */
	{
		uint nrow, ncol;

		count_entries(fname_in, &nrow, &ncol, NULL);

		if(ncol == 0)
		{
			ERR_INTRO "'count_entries' failed\n", ERR_ARG);
			return EXIT_FAILURE;
		}
		double *input = malloc_double(nrow * ncol);

		fpfile2double(fname_in, input, nrow, ncol, 0);

		const uint len = nrow;
		uint tau[ncol >> 1], i;
		double *(alvar[ncol >> 1]),
		       *(bins[ncol >> 1]);

		/**> Main loop starts here */
		show_progress(0, 0,  ncol >> 1);
		for(i = 0; i < ncol; i++)
		{
			if(!(i & 1)) /**< Only odd columns contain data, the even ones hold time stamps */
				continue;

            double *in_cp = malloc_double(len);
			uint j;

			for(j = 0; j < len; j++)
				in_cp[j] = input[j * ncol + i]; /**< This isn't a continuous array, so no 'memcpy' */

			if(use_timeseries_def)
			{
				const double tau_0 = fabs(input[ncol + i - 1] - input[i - 1]);
				alvar[i >> 1] = get_alvar_timeseries(in_cp, len, &tau[i >> 1], tau_0);
			}
			else
				alvar[i >> 1] = get_alvar(in_cp, len, &tau[i >> 1]);

			free(in_cp);
			show_progress(1, 1, 0);

			const double tot_time = fabs(input[(len - 1) * ncol + i - 1] - input[i - 1]);
			bins[i >> 1] = malloc_double(tau[i >> 1]);
			double tt;
			for(j = 0, tt = tot_time / tau[i >> 1]; j < tau[i >> 1]; j++)
				bins[i >> 1][j] = tt * (j + 1.);
		}
		show_progress(0, 2, 0);
		/**< Main loop ends */

		for(i = 1; i < ncol >> 1; i++)
			if(tau[i - 1] != tau[i])
				ERR_INTRO "length of the resulting arrays unequal\n", ERR_ARG);

		free(input);

		double *alvar_s = malloc_double(tau[0] * ncol);

		uint j;
		for(i = 0, j = 0; i < ncol >> 1; i++)
		{
			memcpy(alvar_s + (j++) * tau[0], bins[i], tau[0] * sizeof(double));
			memcpy(alvar_s + (j++) * tau[0], alvar[i], tau[0] * sizeof(double));
			free(bins[i]);
			free(alvar[i]);
		}

		double2fpfile(fname_out, alvar_s, tau[0], ncol);

		free(alvar_s);

		plotdata(fname_out, ncol);
	}
	return EXIT_SUCCESS;
}


uchar parse_args(const int argc, char **argv, char *__restrict fname_in, char *__restrict fname_out,
				const uint len_fname, uchar *__restrict time_data)
{
	int i;
	fprintf(stdout, "alvar read:\n");
	for(i = 0; i < argc; i++)
		fprintf(stdout, "  argv[%i] = '%s'\n", i, argv[i]);

	if(argc != 3 && argc != 4)
	{
		fprintf(stdout,
				"\nHere's the usage:\n" \
				"> alvar TIME_DAT[binary] FILE_IN[string] FILE_OUT[string]\n" \
				"  TIME_DAT tells that the input is a time series.\n" \
				"           Setting it to 0 will change to the original\n" \
				"           definition of the Allan variance of\n" \
				"           sigma^2_y(tau) = .5 * \n" \
				"           <( \\bar y_{n+1} - \\bar y_{n} )^2 )>\n" \
				"           used for frequency data.\n" \
				"  FILE_IN  is the name of the file to be read from.\n" \
				"           Setting it to NULL can be used for\n" \
				"           debugging. The format has to be\n" \
				"           #TIME1 DATA1 TIME2 DATA2 ...\n" \
				"           .1 -1.2984 .4 2.498 ...\n" \
				"           ...\n" \
				"  FILE_OUT is the name of the file that is created\n" \
				"           to store the result in.\n" \
				"           Not providing it will make the computer\n" \
				"           choose a name.\n" \
				"\nExample\n-------\n" \
				"> alvar 1 input.dat output.dat\n" \
				"\nGood luck.\n\n");
		return 0;
	}

	unsigned long int uli = strtoul(argv[1], NULL, 10);
	if(uli != 0 && uli != 1)
	{
		ERR_INTRO "the time/frequency specifier is wrong\n", ERR_ARG);
		return 0;
	}
	else
	{
		*time_data = (uchar)uli;
		fprintf(stdout, "mode: %s\n", *time_data == 1 ? "time" : "frequency");
	}

	strncpy(fname_in, argv[2], len_fname);
	if(!strncmp(fname_in, "NULL", len_fname))
		fprintf(stdout, "debugging mode\n");
	else
	{
		fprintf(stdout, "input: taking '%s'\n", fname_in);
		if(!isfile(fname_in))
		{
			ERR_INTRO "the file is inexistent\n", ERR_ARG);
			return 0;
		}
		if(argc == 4)
		{
			strncpy(fname_out, argv[3], len_fname);
			fprintf(stdout, "output: taking '%s'\n", fname_out);
		}
		else
		{
			strncpy(fname_out, fname_in, len_fname);
			char *pch = strrchr(fname_out, '.');
			const uint pos = pch - fname_out;
			assert('.' == fname_out[pos]);
			snprintf(fname_out + pos, len_fname - pos, "_out.dat");
			fprintf(stdout, "output: taking '%s'\n", fname_out);
		}
	}
	return 1;
}


double get_sum_sqdev(const double *__restrict x, const uint tau, const uint len)
{
	const uint ssq_len = len - tau * 2;
	double res = 0.;
	uint j;
	for(j = 0; j < ssq_len; j++)
	{
		double ssq = 0.;
		uint i;
		for(i = j; i < j + tau; i++)
			ssq += x[i] - x[i + tau];
		ssq /= (double)tau;
		ssq *= ssq;
		res += ssq;
	}
	res /= (double)ssq_len;

	return res;
}


/** \brief Estimates the Allan variance.
 *
 * \param x const double *__restrict The input data.
 * \param len const uint The length of the input data.
 * \param len_res uint *__restrict
 * \return double *
 *
 * This implementation follows the original definition
 * sigma^2_y(tau) = .5 * <( \bar y_{n+1} - \bar y_{n} )^2 )>
 * where
 * \bar y_{n} = 1 / tau \int_{t}^{t + tau} y(t') dt'
 * is the time average.
 */
double *get_alvar(const double *__restrict x, const uint len, uint *__restrict len_res)
{
	const uint len_r = floor(len / 2.) - 1;
	double *res, *x_cp;
	res = malloc_double(len_r);
	x_cp = malloc_double(len);

	memset(res, 0xDE, len_r * sizeof(double));
	memcpy(x_cp, x, len * sizeof(double));

	{
		uint i;
		double mean = 0.;
		for(i = 0; i < len; i++)
			mean += x_cp[i];
		mean /= (double)len;
		for(i = 0; i < len; i++)
			x_cp[i] -= mean;
	}

	uint tau;
	#pragma omp parallel for private(tau), firstprivate(len, len_r), shared(res, x_cp)
	for(tau = 1; tau <= len_r; tau++)
		res[tau - 1] = get_sum_sqdev(x_cp, tau, len) / (2. * tau);

	free(x_cp);
	*len_res = len_r;

	return res;
}


/** \brief Estimates the Allan variance.
 *
 * \param x const double *__restrict The input data.
 * \param len const uint The length of the data.
 * \param len_res uint *__restrict
 * \param tau_0 const double
 * \return double *
 *
 * Estimates the Allan deviation in the formulation
 * for a time series.
 *
 * Ref.: IEEE on Ultrasonics, Ferroelectrics, and Frequency Control,
 * vol. 57, no. 1, Jan 2010, L. Galleani
 */
double *get_alvar_timeseries(const double *__restrict x, const uint len, uint *__restrict len_res, const double tau_0)
{
	uint len_r = floor((len - 1.) / 2.);

	double *res, *x_cp;
	res = malloc_double(len_r);
	x_cp = malloc_double(len);

	memcpy(x_cp, x, len * sizeof(double));

	{
		uint i;
		double mean = 0.;
		for(i = 0; i < len; i++)
			mean += x_cp[i];
		mean /= (double)len;
		for(i = 0; i < len; i++)
			x_cp[i] -= mean;
	}

	uint n;
	#pragma omp parallel for private(n), firstprivate(len, len_r, tau_0), shared(res, x_cp)
	for(n = 1; n < len_r + 1; n++) /* n = 1 ... (N - 1) / 2 */
	{
		uint i;
		double ssq = 0.;
		for(i = 0; i < len - (n << 1); i++)
		{
			double ssq_temp = x_cp[i + (n << 1)] - 2. * x_cp[i + n] + x_cp[i];
			ssq += ssq_temp * ssq_temp;
		}
		double tau = n * tau_0,
		       norm = 2. * tau * tau;
		norm *= (len - 2. * n);
		assert(norm > 0.);

		res[n - 1] = ssq / norm;
	}
	free(x_cp);
	*len_res = len_r;
	return res;
}


uint fpfile2double(const char *__restrict fname, double *__restrict out, const uint row, const uint col, const uchar chatty)
{
	FILE *readfile = fopen(fname, "r");
	if(readfile == NULL)
	{
		ERR_INTRO "'fopen' failed for '%s'\n", ERR_ARG, fname);
		return 0;
	}
	uint i = 0, j = 0, s = 0, ncom = 0;
	char tc;
	const char cmt = '#';
	double temp;
	fpos_t pos;
	fgetpos(readfile, &pos);
	while(fscanf(readfile, "%c", &tc) != EOF)
	{
		if(tc == cmt)
		{
			ncom++;
			while(fscanf(readfile, "%c", &tc) != EOF) /**< Carry on scanning until EOF or line break */
				if(tc == '\n')
					break;
		}
		else if(tc == '0' || tc == '1' || tc == '2' || tc == '3' || tc == '4' ||
				tc == '5' || tc == '6' || tc == '7' || tc == '8' || tc == '9' ||
				tc == '.' || tc == '-' || tc == 'e' || tc == 'E')
		{
			fsetpos(readfile, &pos); /**< Undo the reading in 'while' */
			while(fscanf(readfile, "%18lg", &temp) != EOF)
			{
				out[i * col + j] = temp;
				if(j < (col - 1))
					j++;
				else
				{
					j = 0;
					i++;
				}
				s++;
			}
		}
		fgetpos(readfile, &pos);
	}
	fclose(readfile);
	assert(i == row && s == row * col);
	if(chatty)
		fprintf(stdout,
				"your file '%s'\n"
				"  has %u entries\n"
				"  has %u comment lines\n",
				fname, s, ncom);
	return s;
}


void double2fpfile(const char *__restrict fname, const double *__restrict out, const uint row, const uint col)
{
	FILE *writefile = fopen(fname, "w");
	if(writefile == NULL)
	{
		ERR_INTRO "'fopen' failed for '%s'\n", ERR_ARG, fname);
		return;
	}
	uint i, j;
	const char cmt = '#';

	char tbuf[64];
	{
		time_t rawtime;
		struct tm *timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(tbuf, 64, "%Y%m%d_%H:%M:%S", timeinfo);
	}

	fprintf(writefile, "%c %s {", cmt, tbuf);
	for(i = 0; i < col; i++)
		fprintf(writefile, "%u%c", i, i == col - 1 ? '}' : ' ');
	fprintf(writefile, "\n");

	for(j = 0; j < row; j++)
		for(i = 0; i < col; i++)
			fprintf(writefile, "%g%c",
					out[i * row + j], i == col - 1 ? '\n' : ' ');

	fprintf(writefile, "\n");
	fclose(writefile);
}


/** \brief Counts numeric entries in a file.
 *
 * \param fname char *__restrict The filename.
 * \param nrow uint *__restrict The number of rows of the file.
 * \param ncol uint *__restrict The number of columns.
 * \param cmt_lines uint *__restrict An array of pointers to comment lines.
 * \return void
 *
 * The file has to resemble a matrix. No additional text entries are allowed.
 * The file is checked for consistence during evaluation. The pointer to the current
 * position in the file will be reset to the initial position.
 */
void count_entries(char *__restrict fname,
				uint *__restrict nrow, uint *__restrict ncol,
				uint *__restrict cmt_lines)
{
	uint i = 0, sum = 0, ncom = 0;
	*ncol = 0;
	*nrow = 0; /**< Use this as an error indicator */
	char tc, lbreak = 0, wspace = 0,
	     cmt = '#'; /**< This is the indicator for a comment line in the file */

	FILE *file = fopen(fname, "r");
	if(file == NULL)
	{
		ERR_INTRO "'fopen' failed for '%s'\n", ERR_ARG, fname);
		return;
	};
	while(fscanf(file, "%c", &tc) != EOF)
		if(tc == cmt)
		{
			if(cmt_lines != NULL)
				cmt_lines[ncom] = *ncol; /**< This can be used as a locater for comment lines */
			ncom++;
			while(fscanf(file, "%c", &tc) != EOF) /**< Carry on scanning until EOF or line break */
				if(tc == '\n')
					break;
		}
		else if(tc == '\n')
		{
			if(lbreak)
				break; /**< This happens in case of an additional line break at the end of the file */
			i++; /**< The counter */
			if(!sum)
				*ncol = sum = i; /**< First line break */
			else if((sum += *ncol) != i)
			{
				ERR_INTRO "lines are not of constant length at pos: %u\n", ERR_ARG, i);
				fclose(file);
				return;
			}
			lbreak = 1;
		}
		else if(tc == '0' || tc == '1' || tc == '2' || tc == '3' || tc == '4' ||
				tc == '5' || tc == '6' || tc == '7' || tc == '8' || tc == '9' ||
				tc == '.' || tc == '-' || tc == 'e' || tc == 'E')
			wspace = lbreak = 0; /**< Got an entry */
		else if(tc == ' ' || tc == '\t')
		{
			if(wspace)
			{
				ERR_INTRO "too much whitespace between numbers at pos: %u\n", ERR_ARG, i);
				fclose(file);
				return;
			}
			i++; /**< The counter */
			wspace = 1;
		}
		else if(tc == 13)
			continue; /**< Carriage return on linux */
		else
		{
			ERR_INTRO "last read character '%c' ASCII: %i at pos %u\n", ERR_ARG, tc, tc, i);
			fclose(file);
			return;
		}

	if(!lbreak)
		sum += *ncol; /**< Adds the last column to the total number in case of no final line break */
	fclose(file);
	*nrow = sum / (*ncol);
	fprintf(stdout,
			"your file '%s'\n"
			"  has %u entries\n"
			"  has %u comment lines\n"
			"  form a %u x %u matrix (row x col)\n",
			fname, sum, ncom, *nrow, *ncol);
}


/** \brief Checks whether a file is a file
 *
 * \param szfile const char * The name of the file.
 * \return uchar Returns 1 if the file is a file.
 *
 */
uchar isfile(const char *szfile)
{
	struct stat statBuffer;
	return (stat(szfile, &statBuffer) >= 0 && /**< Make sure it exists */
			!(statBuffer.st_mode & S_IFDIR)); /**< and it's not a directory */
}


/** \brief Prints the progress of the calculation as a bar.
 *
 * \param inc const uint The number of successful runs between two calls of this function.
 * \param order const uchar The order to obey.
 * \param ntot const uint The total number of runs.
 * \return void
 *
 * Initialize with order = 0, update with 1, finish with 2.
 * In accordance with OpenMP specifications.
 */
void show_progressbar(const uint inc, const uchar order, const uint ntot)
{
	static char bar[100], printbar[100];
	static uint sn = 0, stot;
	static double start_t;
	if(order == 0)
	{
		fprintf(stdout, "processing info:\n  %i processors available to this program\n",
				omp_get_num_procs());
		start_t = omp_get_wtime();
		memset(bar, '-', 100);
		stot = ntot;
		fprintf(stdout, "progress (%u steps):\n|%-100.100s|", stot, printbar);
		fflush(stdout);
	}
	else if(order == 1)
	{
		sn += inc;
		assert(printbar != bar);
		assert(sn <= stot);
		memcpy(printbar, bar, sn * 100 / stot);
		fprintf(stdout, "%c|%-100.100s|", 13, printbar);
		fflush(stdout);
	}
	else if(order == 2)
	{
		const uint fin = sn * 100 / stot;
		assert(fin <= 100);
		memcpy(printbar, bar, fin);
		fprintf(stdout, "%c|%-100.100s|", 13, printbar);
		fprintf(stdout, "\nruntime:\n  %.2g s\n",
				fabs(omp_get_wtime() - start_t));
		sn = 0;
	}
	else
		ERR_INTRO "wrong order code, doing nothing\n", ERR_ARG);
}

/** \brief Prints the progress of the calculation in a percentage.
 *
 * \param inc const uint The number of successful runs between two calls of this function.
 * \param order const uchar The order to obey.
 * \param ntot const uint The total number of runs.
 * \return void
 *
 * Initialize with order = 0, update with 1, finish with 2.
 * In accordance with OpenMP specifications.
 */
void show_progress(const uint inc, const uchar order, const uint ntot)
{
	static uint sn = 0, stot;
	static double start_t;
	if(order == 0)
	{
		fprintf(stdout, "processing info:\n  %i processors available to this program\n",
				omp_get_num_procs());
		start_t = omp_get_wtime();
		stot = ntot;
		fprintf(stdout, "progress (%u steps):\n  %3u %%", stot, 0u);
		fflush(stdout);
	}
	else if(order == 1)
	{
		sn += inc;
		fprintf(stdout, "%c  %3u %%", 13, sn * 100 / stot);
		fflush(stdout);
	}
	else if(order == 2)
	{
		const uint fin = sn * 100 / stot;
		assert(fin <= 100);
		fprintf(stdout, "%c  %3u %%", 13, fin);
		fprintf(stdout, "\nruntime:\n  %.2g s\n",
				fabs(omp_get_wtime() - start_t));
		sn = 0;
	}
	else
		ERR_INTRO "wrong order code, doing nothing\n", ERR_ARG);
}


/** \brief Allocates and initializes an array of memory.
 *
 * \param size const uint The length of the array.
 * \return double * The allocated and initialized array.
 *
 */
double *malloc_double(const uint size)
{
	double *arr = malloc(size * sizeof(double));
	if(arr == NULL)
	{
		ERR_INTRO "'malloc' failed\n", ERR_ARG);
		exit(EXIT_FAILURE);
	}
	memset(arr, 0xFF, size * sizeof(double));
	return arr;
}


char *replace_char(const char *str_in, const char find, const char replace)
{
	char *ret = strdup(str_in),
	     *s = strdup(str_in);
	while(*s != 0)
	{
		if(*s == find)
		{
			*ret++ = replace;
			++s;
		}
		else
			*ret++ = *s++;
	}
	*ret = '\0';
	return ret;
}


void plotdata(const char *__restrict fname, const uint ncol)
{
	const char *cmdtmp = "this_is_just_a_temporary_thing";

	FILE *gnufile = fopen(cmdtmp, "w");
	if(gnufile == NULL)
	{
		ERR_INTRO "can't open file", ERR_ARG);
		exit(EXIT_FAILURE);
	}

	char timebuf[20], *plttitle;
	{
		time_t tt;
		time(&tt);
		strftime(timebuf, 20, "%Y-%m-%d_%H-%M-%S", localtime(&tt));
		plttitle = replace_char(timebuf, '_', ' ');
	}

	fprintf(gnufile,
			"set term postscript enhanced color dl 1\n" \
			"set key top Left nobox\n" \
			"set grid back\n" \
			"set xr [:]\n" \
			"set yr [:]\n" \
			"set log xy\n" \
			"set format y '%%g'\n" \
			"set format x '%%g'\n" \
			"set ylabel '{/Symbol s}_y^2'\n" \
			"set xlabel 'Bin length'\n" \
			"set title '%s'\n" \
			"set out '%s.ps'\n" \
			"plot ", plttitle, timebuf);

	uint i;
	for(i = 0; i < ncol >> 1; i++)
		fprintf(gnufile, "'%s' every ::2 using %u:%u w l t 'Column %u'%c",
				fname, (i << 1) + 1, (i << 1) + 2, i + 1, i == (ncol >> 1) - 1 ? '\n' : ',');

	fclose(gnufile);

	fprintf(stdout, "gnuplot message: ");
	fflush(stdout);
	i = system("gnuplot this_is_just_a_temporary_thing");
	if(i)
		fprintf(stderr, "\ngnuplot command returned '1' -- the plotting went wrong!\n");
	else
		fprintf(stdout, "plotted '%s.ps'\n", timebuf);

	remove(cmdtmp);
}
