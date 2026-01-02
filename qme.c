#include "qme.h"

int main(int argc, char* argv[])
{
	t_qme *qme = malloc(sizeof(t_qme));
	FILE *opt_file = NULL;

	opt_file = sfopen(argv[1], 0);
	get_option(qme, opt_file);

#if 1
	switch(qme->jtype){
		case 0: fret_prep(qme);
				dyn_inco(qme);
				break;
		case 1: mrt_prep(qme);
				dyn_inco(qme);
				break;
		case 10: fret_prep(qme);
				dyn_inco(qme);
				break;
		case 11: mrt_prep(qme);
				dyn_inco(qme);
				break;
		case 20: fret_prep(qme);
				comp_tr(qme);
				break;
		case 21: mrt_prep(qme);
				comp_tr(qme);
				break;
		case 30: avg_split(qme);
		default: break;
	}
#endif

	fclose(opt_file);
	return 0;
}
