#include "util.h"

FILE* sfopen(char *filename, int mode)
{
	FILE *checkfile;

	if(mode == 0) checkfile = fopen(filename, "r");
	else if(mode == 1) checkfile = fopen(filename, "w");
	else if(mode == 2) checkfile = fopen(filename, "a");
	else if(mode == 3) checkfile = fopen(filename, "r+");
	else if(mode == 4) checkfile = fopen(filename, "w+");
	else {
		printf("ERROR:: unknown error in file opening\n");
		exit(-1);
	}

	if(checkfile == NULL) {
		printf("ERROR::%s already exists or cannot be opened\n", filename);
		exit(-1);
	}

	return(checkfile);
}

void *save_calloc(char *name, int nelem, double *mem, size_t size)
{
	void *p;                                                                                          
	p = NULL;

	if((nelem == 0) || (size == 0)) p = NULL;                                      
	else p = calloc((size_t)nelem, size);
	*mem += (double)(nelem * (int)size);

	return p; 
}

void *save_malloc(char *name, int nelem, double *mem, size_t size)
{
	void *p;
	long long int size_tot = (long long int)nelem * (long long int)size;
	p = NULL;

	if(size_tot == 0) p = NULL;
	else {
		if((p = malloc(size_tot))==NULL) {
			fprintf(stderr, "Not enough memory. Failed to malloc %d bytes for %s\n", size_tot, name);
			exit(EXIT_FAILURE);
		}
		(void) memset (p, 0, size_tot);
	}
	*mem += (double)(nelem * (int)size);

	return p;
}

int get_lines(char *filename)
{
	int i;
	int lines = 0;
	char buf[255];

	FILE *f = sfopen(filename, 0);
	while(fgets(buf, 255, f) != NULL) lines++;
	fclose(f);

	return lines;
}

void upper2low(char *line)
{
	int i, len;
	len = strlen(line);
	for (i=0; i<len; i++) line[i] = (char)tolower(line[i]);
}

void prt_vib_dmat(int ndim, int nvib, double *vib_mat, double fac)
{
	int i, j;

	for(i=0; i<ndim; i++){
		if(i % nvib == 0){
			for(j=0; j<ndim; j++){
				if(j % nvib == 0) fprintf(stderr, "|");
				fprintf(stderr, "---------------");
			}
			fprintf(stderr, "|\n");
		}
		for(j=0; j<ndim; j++){
			if(j % nvib == 0) fprintf(stderr, "|");
			if(i != j && fabs(vib_mat[i+ndim*j]) < SMALL) fprintf(stderr, "               ");
			else fprintf(stderr, "%15.8lf", vib_mat[i+ndim*j] * fac);
		}
		fprintf(stderr, "|\n");
	}
	for(i=0; i<ndim; i++){
		if(i % nvib == 0) fprintf(stderr, "|");
		fprintf(stderr, "---------------");
	}
	fprintf(stderr, "|\n");
}

void prt_vib_dmat_sci(int ndim, int nvib, double *vib_mat, double fac)
{
	int i, j;

	for(i=0; i<ndim; i++){
		if(i % nvib == 0){
			for(j=0; j<ndim; j++){
				if(j % nvib == 0) fprintf(stderr, "|");
				fprintf(stderr, "---------------");
			}
			fprintf(stderr, "|\n");
		}
		for(j=0; j<ndim; j++){
			if(j % nvib == 0) fprintf(stderr, "|");
			fprintf(stderr, "%15.8le", vib_mat[i+ndim*j] * fac);
		}
		fprintf(stderr, "|\n");
	}
	for(i=0; i<ndim; i++){
		if(i % nvib == 0) fprintf(stderr, "|");
		fprintf(stderr, "---------------");
	}
	fprintf(stderr, "|\n");
}

void prt_vib_zmat(int ndim, int nvib, complex *vib_mat, double fac)
{
	int i, j;

	for(i=0; i<ndim; i++){
		if(i % nvib == 0){
			for(j=0; j<ndim; j++){
				if(j % nvib == 0) fprintf(stderr, "|");
				fprintf(stderr, "---------------");
			}
			fprintf(stderr, "|\n");
		}
		for(j=0; j<ndim; j++){
			if(j % nvib == 0) fprintf(stderr, "|");
			if(i != j && fabs(vib_mat[i+ndim*j].real) < SMALL) fprintf(stderr, "               ");
			else fprintf(stderr, "%15.8lf", vib_mat[i+ndim*j].real * fac);
		}
		fprintf(stderr, "|\n");
	}
	for(i=0; i<ndim; i++){
		if(i % nvib == 0) fprintf(stderr, "|");
		fprintf(stderr, "---------------");
	}
	fprintf(stderr, "|\n");
}

void prt_vib_zmat_sci(complex *vib_mat, int ndim, int nvib, double fac)
{
	int i, j;

	for(i=0; i<ndim; i++){
		if(i % nvib == 0){
			for(j=0; j<ndim; j++){
				if(j % nvib == 0) fprintf(stderr, "|");
				fprintf(stderr, "---------------");
			}
			fprintf(stderr, "|\n");
		}
		for(j=0; j<ndim; j++){
			if(j % nvib == 0) fprintf(stderr, "|");
			if(i != j && fabs(vib_mat[i+ndim*j].real) < SMALL) fprintf(stderr, "               ");
			else fprintf(stderr, "%15.8le", vib_mat[i+ndim*j].real * fac);
		}
		fprintf(stderr, "|\n");
	}
	for(i=0; i<ndim; i++){
		if(i % nvib == 0) fprintf(stderr, "|");
		fprintf(stderr, "---------------");
	}
	fprintf(stderr, "|\n");
}

void prt_vib_zmat_imag(int ndim, int nvib, complex *vib_mat, double fac)
{
	int i, j;

	for(i=0; i<ndim; i++){
		if(i % nvib == 0){
			for(j=0; j<ndim; j++){
				if(j % nvib == 0) fprintf(stderr, "|");
				fprintf(stderr, "---------------");
			}
			fprintf(stderr, "|\n");
		}
		for(j=0; j<ndim; j++){
			if(j % nvib == 0) fprintf(stderr, "|");
			if(i != j && fabs(vib_mat[i+ndim*j].imag) < SMALL) fprintf(stderr, "               ");
			else fprintf(stderr, "%15.8lf", vib_mat[i+ndim*j].imag * fac);
		}
		fprintf(stderr, "|\n");
	}
	for(i=0; i<ndim; i++){
		if(i % nvib == 0) fprintf(stderr, "|");
		fprintf(stderr, "---------------");
	}
	fprintf(stderr, "|\n");
}
