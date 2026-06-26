#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#define MAXBUF 3000
#define MAXDIG 30

int main(int argc, char *argv[])
{
	FILE *fin1;
	char buffer[MAXBUF];
	char judge;
	char sitenu[MAXDIG], synsit[MAXDIG], synchg[MAXDIG], synsub[MAXDIG];
	char nonsit[MAXDIG], nonchg[MAXDIG], nonsub[MAXDIG];
	int aaa, aab, aac, aad, siten, number;
	double synsite, synchge, synsubs, nonsite, nonchge, nonsubs;
	double syncase, noncase, tsynsub, tnonsub, tsynsqu, tnonsqu, syndisc, nondisc;
	double purify, nosele, advant, juddisc, discoun, tvalue;
	double avsynsb, avnonsb, avsynsq, avnonsq, sesynsb, senonsb, exceede, stopcod, recalcu;
	char flgints, flgintn, flgsel;
	double zza, zzb, zzc, intschg, intnchg, totint, prob, negprob, posprob, totprob, noselect, negasele, posisele;
	double expschg, expnchg, probsyn, probnon, maxprob;
	double logprob;
	int equal, dsdn, dnds, positive, negative, exclude;





/****** one-tail negative ******/
	int neg1tail;


	neg1tail=0;
/****** one-tail negative ******/







	if(argc!=2) /* Usage */
	{
		printf("adaptsite-t : a program to test neutrality at single amino acid sites\n"); fflush(stdout);
		printf("Use         : adaptsite-t filename\n"); fflush(stdout);
		printf("filename    : output file from adaptsite-p or adaptsite-d\n"); fflush(stdout);
		exit(1);
	}
	if((fin1=fopen(argv[1], "r"))==NULL) /* Filename */
	{
		printf("File could not be found\n"); fflush(stdout);
		exit(1);
	}
/******
	number=atoi(argv[2]);
******/

	recalcu=0.0; /****** The number of re-calculation ******/
	syncase=noncase=0.0; /****** The number of codon sites where synonymous or nonsynonymous site exists ******/
	tsynsub=tnonsub=0.0; /****** The total number of synonymous or nonsynonymous substitutions ******/
	tsynsqu=tnonsqu=0.0; /****** The total number of square of synonymous or nonsynonymous substitutions ******/
	syndisc=nondisc=0.0; /****** The number of codon sites where synonymous or nonsynonymous site was absent ******/
	purify=nosele=advant=0.0; /****** The number of codon sites where purifying, no, or advantageous selection is operating ******/
	juddisc=0.0; /****** The number of sites where judge was not done ******/
	discoun=0.0; /****** The number of sites where any number was not calculated due to following two reasons ******/
	exceede=0.0; /****** The number of sites where the number of possible assignments exceeded 10,000 ******/
	stopcod=0.0; /****** The number of sites where the assignment was not available due to the stop codon ******/
	synsite=nonsite=0.0; /****** The number of synonymous or nonsynonymous sites ******/
	synchge=synchge=0.0; /****** The number of synonymous or nonsynonymous change ******/
	synsubs=nonsubs=0.0; /****** The number of synonymous or nonsynonymous substitutions ******/
	avsynsb=avnonsb=0.0; /****** The number of average synonymous or nonsynonymous substitutions ******/
	avsynsq=avnonsq=0.0; /****** The number of average square of synonymous or nonsynonymous substitutions ******/
	sesynsb=senonsb=0.0; /****** The number of standard variation of synonymous or nonsynonymous substitutions ******/
	tvalue=0.0; /****** t-value for judgement ******/
	aad=0; /****** The number of iterations ******/
	equal=dsdn=dnds=positive=negative=exclude=0;






/****** one-tailed negative ******/
	printf("Cod.#\tone-tail-posp\tone-tail-negp\tone-tail-nsel\ttwo-tail-p\tSelect\ttwo-tail-1-p\n"); fflush(stdout);
/****** one-tailed negative ******/







	while(fgets(buffer, MAXBUF, fin1)!=NULL)
	{
		if(isdigit(buffer[0])) /* Site number */
		{
			aaa=0;
			aab=0;
			for(aaa; !isspace(buffer[aaa]); aaa++)
			{
				sitenu[aab]=buffer[aaa];
				++aab;
			}
			sitenu[aab]='\0';
			siten=atoi(sitenu);
			printf("%d\t", siten); fflush(stdout);

			for(aaa; isspace(buffer[aaa]); aaa++) { ; } /****** skipping spaces ******/
			if(isdigit(buffer[aaa])) /****** if results exist ******/
			{
				for(aaa; !isspace(buffer[aaa]); aaa++) { ; } /****** skipping assignment number ******/
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /******* skipping spaces ******/

/******/
				if(buffer[aaa]=='n')
				{
					printf("excluded\n"); fflush(stdout);
					++exclude;
					goto tsugi2;
				}
/******/

				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /****** synonymous sites ******/
				{
					synsit[aab]=buffer[aaa];
					++aab;
				}
				synsit[aab]='\0';
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /* Skipping the spaces */
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /* synonymous changes */
				{
					synchg[aab]=buffer[aaa];
					++aab;
				}
				synchg[aab]='\0';
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /* Skipping the spaces */
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /* synonymous substitutions per site */
				{
					synsub[aab]=buffer[aaa];
					++aab;
				}
				synsub[aab]='\0';
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /* Skipping the spaces */
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /****** nonsynonymous site ******/
				{
					nonsit[aab]=buffer[aaa];
					++aab;
				}
				nonsit[aab]='\0';
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /* Skipping the spaces */
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /* nonsynonymous changes */
				{
					nonchg[aab]=buffer[aaa];
					++aab;
				}
				nonchg[aab]='\0';
				for(aaa; isspace(buffer[aaa]); aaa++) { ; } /* Skipping the spaces */
				aab=0;
				for(aaa; !isspace(buffer[aaa]); aaa++) /* nonsynonymous substitutions per site */
				{
					nonsub[aab]=buffer[aaa];
					++aab;
				}
				nonsub[aab]='\0';

				synsite=atof(synsit);
				synchge=atof(synchg);
				synsubs=atof(synsub); /* Convert number of synonymous substitutions */
				nonsite=atof(nonsit);
				nonchge=atof(nonchg);
				nonsubs=atof(nonsub); /* Convert number of nonsynonymous substitutions */

				if((synsubs!=9999.0 && nonsubs!=9999.0) && (!(synsubs==0.0 && nonsubs==0.0))) /* Count judgement */
				{
/******
					tvalue=((synsubs-nonsubs)/(sqrt((synsubs/synsite)+(nonsubs/nonsite))));
					if(tvalue>1.96)
					{
						printf("\t\t%.10f \tP\t", tvalue); fflush(stdout);
					}
					else if((-1.96<=tvalue) && (tvalue<=1.96))
					{
						printf("\t\t%.10f \tN\t", tvalue); fflush(stdout);
					}
					else if(tvalue<-1.96)
					{
						printf("\t\t%.10f \tA\t", tvalue); fflush(stdout);
					}
******/
/****** ******/ 
flgints='n'; /****** converting raw number of synonymous changes to integer ******/
for(zza=0.0; zza<=10000.0; zza=zza+1.0)
{
	if(zza-0.5<=synchge && synchge<zza+0.5)
	{
		intschg=zza;
		flgints='y';
	}
}
if(flgints=='n')
{
	printf("ERROR 1\n"); fflush(stdout);
	exit(1);
}

flgintn='n'; /****** converting raw number of nonsynonymous changes to integer ******/
for(zza=0.0; zza<=10000.0; zza=zza+1.0)
{
	if(zza-0.5<=nonchge && nonchge<zza+0.5)
	{
		intnchg=zza;
		flgintn='y';
	}
}
if(flgintn=='n')
{
	printf("ERROR 2\n"); fflush(stdout);
	exit(1);
}
totint=intschg+intnchg;
expschg=totint*(synsite/3.0);
expnchg=totint*(nonsite/3.0);
probsyn=synsite/3.0;
probnon=nonsite/3.0;





/******/
if(totint==0.0)
{
	goto koko1;
}
/******/






/****** One-tailed ******/
					negprob=0.0;
					posprob=0.0;
					for(zza=intschg; zza<=totint; zza=zza+1.0) /****** First tail ******/
					{
						zzc=totint;
						prob=1.0;
/******/
						logprob=0.0;
/******/
						for(zzb=zza; zzb>=1.0; zzb=zzb-1.0) 
						{
/******
							prob=prob*(zzc/zzb); 
******/
/******/
							logprob=logprob+log10(zzc)-log10(zzb);
/******/
							zzc=zzc-1.0;
						}
/******
						prob=prob*pow(probsyn, zza)*pow(probnon, totint-zza);
******/
/******/
						logprob=logprob+log10(probsyn)*zza+log10(probnon)*(totint-zza);
						prob=pow(10.0, logprob);
/******/
						if(zza==intschg)
						{
							maxprob=prob;
						}
						negprob=negprob+prob;
					}
					for(zza=intnchg; zza<=totint; zza=zza+1.0)
					{
						zzc=totint;
						prob=1.0;
/******/
						logprob=0.0;
/******/
						for(zzb=zza; zzb>=1.0; zzb=zzb-1.0)
						{
/******
							prob=prob*(zzc/zzb);
******/
/******/
							logprob=logprob+log10(zzc)-log10(zzb);
/******/
							zzc=zzc-1.0;
						}
/******
						prob=prob*pow(probnon, zza)*pow(probsyn, totint-zza);
******/
/******/
						logprob=logprob+log10(probnon)*zza+log10(probsyn)*(totint-zza);
						prob=pow(10.0, logprob);
/******/
						if(zza==intnchg)
						{
							maxprob=prob;
						}
						posprob=posprob+prob;
					}






					printf("%.10f\t%.10f\t", posprob, negprob); fflush(stdout);





/****** one-tailed negative ******/
					if(negprob<0.05)
					{
						printf("P\t"); fflush(stdout);
						++neg1tail;
					}
					else if(negprob>=0.05)
					{
						printf("N\t"); fflush(stdout);
					}
/****** one-tailed negative ******/




/****** Two-tailed ******/
					if(intschg==expschg) /****** As expected ******/
					{
						flgsel='e';
						totprob=1.0;
					}
					else if(intschg>expschg) /****** Negative selection ******/
					{
						flgsel='n';
						totprob=0.0;
						for(zza=intschg; zza<=totint; zza=zza+1.0) /****** First tail ******/
						{
							zzc=totint;
							prob=1.0;
/******/
							logprob=0.0;
/******/
							for(zzb=zza; zzb>=1.0; zzb=zzb-1.0) 
							{ 
/******
								prob=prob*(zzc/zzb); 
******/
/******/
								logprob=logprob+log10(zzc)-log10(zzb);
/******/
								zzc=zzc-1.0;
							}
/******
							prob=prob*pow(probsyn, zza)*pow(probnon, totint-zza);
******/
/******/
							logprob=logprob+log10(probsyn)*zza+log10(probnon)*(totint-zza);
							prob=pow(10.0, logprob);
/******/
							if(zza==intschg)
							{
								maxprob=prob;
							}
							totprob=totprob+prob;
						}
						for(zza=0.0; zza<intschg; zza=zza+1.0) /****** Second tail ******/
						{
							zzc=totint;
							prob=1.0;
/******/
							logprob=0.0;
/******/
							if(zza==0.0)
							{
/******
								prob=prob*pow(probnon, totint);
******/
/******/
								logprob=logprob+log10(probnon)*totint;
								prob=pow(10.0, logprob);
/******/
							}
							else if(zza!=0.0)
							{
								for(zzb=zza; zzb>=1.0; zzb=zzb-1.0)
								{
/******
									prob=prob*(zzc/zzb);
******/
/******/
									logprob=logprob+log10(zzc)-log10(zzb);
/******/
									zzc=zzc-1.0;
								}
/******
								prob=prob*pow(probsyn, zza)*pow(probnon, totint-zza);
******/
/******/
								logprob=logprob+log10(probsyn)*zza+log10(probnon)*(totint-zza);
								prob=pow(10.0, logprob);
/******/
							}
							if(prob<=maxprob)
							{
								totprob=totprob+prob;
							}
							else if(prob>maxprob)
							{
								;
							}
						}
					}
					else if(intschg<expschg) /****** Positive selection ******/
					{
						flgsel='p';
						totprob=0.0;
						for(zza=intnchg; zza<=totint; zza=zza+1.0)
						{
							zzc=totint;
							prob=1.0;
/******/
							logprob=0.0;
/******/
							for(zzb=zza; zzb>=1.0; zzb=zzb-1.0)
							{
/******
								prob=prob*(zzc/zzb);
******/
/******/
								logprob=logprob+log10(zzc)-log10(zzb);
/******/
								zzc=zzc-1.0;
							}
/******
							prob=prob*pow(probnon, zza)*pow(probsyn, totint-zza);
******/
/******/
							logprob=logprob+log10(probnon)*zza+log10(probsyn)*(totint-zza);
							prob=pow(10.0, logprob);
/******/
							if(zza==intnchg)
							{
								maxprob=prob;
							}
							totprob=totprob+prob;
						}
						for(zza=0.0; zza<intnchg; zza=zza+1.0) /****** Second tail ******/
						{
							zzc=totint;
							prob=1.0;
/******/
							logprob=0.0;
/******/
							if(zza==0.0)
							{
/******
								prob=prob*pow(probsyn, totint);
******/
/******/
								logprob=logprob+log10(probsyn)*totint;
								prob=pow(10.0, logprob);
/******/
							}
							else if(zza!=0.0)
							{
								for(zzb=zza; zzb>=1.0; zzb=zzb-1.0)
								{
/******
									prob=prob*(zzc/zzb);
******/
/******/
									logprob=logprob+log10(zzc)-log10(zzb);
/******/
									zzc=zzc-1.0;
								}
/******
								prob=prob*pow(probnon, zza)*pow(probsyn, totint-zza);
******/
/******/
								logprob=logprob+log10(probnon)*zza+log10(probsyn)*(totint-zza);
								prob=pow(10.0, logprob);
/******/
							}
							if(prob<=maxprob)
							{
								totprob=totprob+prob;
							}
							else if(prob>maxprob)
							{
								;
							}
						}
					}






					if(totprob<0.05)
					{
						if(flgsel=='n')
						{
							printf("%.10f\tP", -1.0*totprob); fflush(stdout); /****** Purifying selection ******/
							printf("\t%.10f\n", -1.0*(1.0-totprob)); fflush(stdout);
							++dsdn;
							++negative;
						}
						else if(flgsel=='p')
						{
							printf("%.10f\tA", totprob); fflush(stdout); /****** Advantageous selection ******/
							printf("\t%.10f\n", 1.0-totprob); fflush(stdout);
							++dnds;
							++positive;
						}
					}
					else if(totprob>=0.05)
					{
						if(flgsel=='n')
						{
							printf("%.10f\tN", -1.0*totprob); fflush(stdout);
							printf("\t%.10f\n", -1.0*(1.0-totprob)); fflush(stdout);
							++dsdn;
						}
						else if(flgsel=='p')
						{
							printf("%.10f\tN", totprob); fflush(stdout);
							printf("\t%.10f\n", 1.0-totprob); fflush(stdout);
							++dnds;
						}
						else if(flgsel=='e')
						{
							printf("+-1.0\t\tN"); fflush(stdout);
							printf("\t0.0\n"); fflush(stdout);
							++equal;
						}
					}






/****** ******/

				}
				else if(synsubs==9999.0 || nonsubs==9999.0)
				{
					++exclude;
					if(synsite==9999.0)
					{
						printf("Gap included\n"); fflush(stdout);
						goto tsugi2;
					}
					printf("Site was absent\n"); fflush(stdout);
				}
				else if(synsubs==0.0 && nonsubs==0.0)
				{






/******/
					koko1: ;
/******/







					++equal;
/******
					printf("t-incalculable\t\t9999.0\tN"); fflush(stdout);
******/
					printf("nosubstitution"); fflush(stdout);
					printf("\n"); fflush(stdout);
				}
			}
			else if(isalpha(buffer[aaa])) /* If results do not exist */
			{
				++exclude;
				for(aaa; buffer[aaa]!='\n'; aaa++)
				{
					if(buffer[aaa]=='0') /****** If at least one 0 existed in the string ******/
					{
						printf("More than 10,000\n"); fflush(stdout);
						goto tsugi2;
					}
					if(buffer[aaa]=='p') /****** If at least one p existed in the string ******/
					{
						printf("Stop codon\n"); fflush(stdout);
						goto tsugi2;
					}
				}
				tsugi2: ;
			}
		}
	}
	fclose(fin1);

	printf("\n"); fflush(stdout);
	printf("Equal      %d\n", equal); fflush(stdout);
	printf("dN>dS      %d\n", dnds); fflush(stdout);
	printf("Positive   %d\n", positive); fflush(stdout);
	printf("dN<dS      %d\n", dsdn); fflush(stdout);
	printf("Negative   %d\n", negative); fflush(stdout);
	printf("Neg1tail   %d\n", neg1tail); fflush(stdout);
	printf("Excluded   %d\n", exclude); fflush(stdout);
	printf("Total      %d\n", equal+dnds+dsdn+exclude); fflush(stdout);
}
