struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct devent *nextde;
} ;
struct c_params {
	int npop;
	int nsam;
	int *config; /* lh3? what is this */
	double **mig_mat;
	double r; /* lh3: recombination rate */
	int nsites; /* lh3: number of sites (for -r option) */
	int is_lh3; /* lh3: added by me */
	double f;
	double track_len;
	double *size; /* lh3: (constant) size of population in history */
	double *alphag; /* lh3: exponential parameter, which is combined with *size */
	struct devent *deventlist ;
} ;
struct m_params {
	double theta;
	int segsitesin;
	int treeflag;
	int treelengthflag;
	int mfreq;
} ;
struct params { 
	struct c_params cp;
	struct m_params mp;
	int commandlineseedflag ;
};
