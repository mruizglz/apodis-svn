/*

Program implemented by Juan Manuel Lopez, Universidad Politecnica de Madrid
july-Augoust 2011
Implementation of Apodis 2

Functions for load setup values


*/

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/resource.h>

#define DEFAULT_CONFIG_FILE "./apodis.conf"

#ifdef DEBUG
#define PDEBUG(x...) do {fprintf(stdout,x);}while(0)
#else
#define PDEBUG(x...)
#endif

char     *conf_file_name;
double   conf_Threshold = -1100000;
double   conf_Sampling = 0.001;

int      conf_Npoints = 32;
int      conf_Nsignals= 7;
int      conf_NcalculateSignals = 2;
int      conf_NModels = 3;

char     *conf_PathMax = "./trainingXX/MaxCoefficients.txt";
char     *conf_PathMin = "./trainingXX/MinCoefficients.txt";

char     *conf_PathModel  =  "./trainingXX/ModelDescription.txt";

//char     *conf_PathR   =  "./trainingXX/RCoefficients.txt";
char     *conf_PathR   =  "./trainingXX/M.txt";

char     *conf_PathNormalize   =  "./trainingXX/Normalize.txt";





static void c_set_string(char *v1, const char *v2, void *t);
static void c_set_int(char *v1, const char *v2, void *t);
static void c_set_unity(char *v1, const char *v2, void *t);
static void c_set_double(char *v1, const char *v2, void *t);


struct ccommand {
	const char *name;
	const int type;
	void (*action) (char *, const char *, void *);
	void *object;
};

typedef struct ccommand Command;


/* Help keep the table below compact */
#define STMT_NO_ARGS 1
#define STMT_ONE_ARG 2

#define S0A STMT_NO_ARGS
#define S1A STMT_ONE_ARG

struct ccommand clist[] = {
	{"Threshold", 		S1A,c_set_double,	&conf_Threshold},
	{"Maximum", 	        S1A,c_set_string,	&conf_PathMax},
	{"Minimum", 	        S1A,c_set_string,	&conf_PathMin},
	{"PathModel", 		S1A,c_set_string,	&conf_PathModel},
	{"PathR", 		S1A,c_set_string,	&conf_PathR},
	{"PathNormalize", 		S1A,c_set_string,	&conf_PathNormalize},
	{"Sampling", 		S1A,c_set_double,	&conf_Sampling},
	{"Npoints",	        S1A,c_set_int,		&conf_Npoints},
	{"Nsignals",	        S1A,c_set_int,		&conf_Nsignals},
	{"NcalculateSignals",	S1A,c_set_int,		&conf_NcalculateSignals},
	{"NModels",	        S1A,c_set_int,		&conf_NModels},
};

static void c_set_string(char *v1, const char *v2, void *t)
{
	if (t)
	{
		PDEBUG("STRING TO STORE: %s\n", v1);
		if (*v1=='"')
		{
			v1++;
			if (v1[strlen(v1)-1]!='"')
			{
				PDEBUG("Fail to close '\"'\n");
				exit(EXIT_FAILURE);
			}
			v1[strlen(v1)-1]='\0';
		}
		//if (*(char **) t != NULL)
		//   free(*(char **) t);
		*(char **) t = (char *)strdup(v1);
		if (!*(char **) t)
		{
			PDEBUG("Unable to strdup in c_set_string\n");
			exit(EXIT_FAILURE);
		}
		//done
		PDEBUG("Config string: %s=%s\n", v2, *(char**)t);
	}
	else
	{
		//skipped
		//PDEBUG("SKIPPED STRING: %s\n", *(char **)t);
	}
}

static void c_set_int(char *v1, const char *v2, void *t)
{
	char *endptr;
	int i;

	if (t)
	{
		i = strtol(v1, &endptr, 0); /* Automatic base 10/16/8 switching */
		if (*v1 != '\0' && *endptr == '\0')
		{
			*(int *) t = i;
			//converted
			PDEBUG("Config int: %s=%d\n", v2, *(int*)t);
		}
		else
		{
			/* XXX should tell line number to user */
			PDEBUG("Unable to convert int in c_set_int\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//skipped
		//PDEBUG("SKIPPED INT: %d\n", *(int*)t);
	}
}

static void c_set_double(char *v1, const char *v2, void *t)
{
	char *endptr;
	double i;

	if (t)
	{
	        i = strtod (v1, &endptr);
		if (*v1 != '\0' && *endptr == '\0')
		{
			*(double *) t = i;
			//converted
			PDEBUG("Config int: %s=%f\n", v2, *(double*)t);

		}
		else
		{
			/* XXX should tell line number to user */
			PDEBUG("Unable to convert int in c_set_int\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		//skipped
		PDEBUG("SKIPPED FLOAT: %d\n", *(int*)t);
	}
}






static void c_set_unity(char *v1, const char *v2, void *t)
{
	if (t)
	{
		PDEBUG("Set %s\n",v2);
		*(int *) t = 1;
	}
}

struct ccommand *lookup_keyword(char *c)
{
	struct ccommand *p;
	for (p = clist;
		p < clist + (sizeof (clist) / sizeof (struct ccommand)); p++)
		{
			if (strcasecmp(c, p->name) == 0)
				return p;
		}
	return NULL;
}

static void apply_command(Command * p, char *args)
{

	switch (p->type)
	{
	case STMT_NO_ARGS:
		(p->action) (NULL, p->name, p->object);
		break;
	case STMT_ONE_ARG:
		(p->action) (args, p->name, p->object);
		break;
	default:
		exit(EXIT_FAILURE);
	}
}

static void trim(char *s)
{
	char *c = s + strlen(s) - 1;

	while (isspace(*c) && c > s)
	{
		*c = '\0';
		--c;
	}
}

static void parse(FILE * f)
{
	char buf[1025], *b, *c;
	Command *p;
	int line = 0;

	while (fgets(buf, 1024, f) != NULL)
	{
		++line;
		b=buf;
		while (isspace(*b))
			++b;

		if (b[0] == '\0' || b[0] == '#' || b[0] == '\n')
			continue;
		/* kill the linefeed and any trailing whitespace */
		trim(b);
		if (b[0] == '\0')
			continue;

		/* look for multiple arguments */
		c = b;
		while (!isspace(*c))
			++c;

		if (*c == '\0')
		{
			/* no args */
			c = NULL;
		}
        else
		{
			/* one or more args */

			*c = '\0';
			++c;
		}

		p = lookup_keyword(b);

		while (isspace(*c))
			++c;

		if (!p)
		{
			PDEBUG( "Line %d: Did not find keyword \"%s\"\n", line,b);
			printf( "Line %d: Did not find keyword \"%s\"\n", line,b);
			exit(EXIT_FAILURE);
		}
		else
		{
/*			if (c==NULL)
				PDEBUG("Found keyword %s in \"%s\" (NULL)!\n",
						p->name, buf);
			else
				PDEBUG("Found keyword %s in \"%s\" (%s)!\n",
						p->name, buf, c);*/

			apply_command(p, c);
		}
	}
}

/*
 * Name: read_config_files
 *
 * Description: Reads config files, then makes sure that
 * all required variables were set properly.
 */
void read_config_files(void)
{
	FILE *config;

	printf ("Reading Configuration FILES \n");

	if (!conf_file_name)
	{
		conf_file_name = DEFAULT_CONFIG_FILE;
	}

	config = fopen(conf_file_name, "r");
	if (!config)
	{
	        printf ("Could not open %s for reading.\n", conf_file_name);
		PDEBUG("Could not open %s for reading.\n", conf_file_name);
		exit(EXIT_FAILURE);
	}
	parse(config);
	fclose(config);
	printf("Configuration FILES were read \n");
}
