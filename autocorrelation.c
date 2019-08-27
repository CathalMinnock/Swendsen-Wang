#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv) {
	int option_index = 0;
	char *input_filename, *output_filename;
	int tau;
	while((option_index = getopt(argc, argv, "i:o:t:" )) != -1) {
		switch(option_index) {
			case 'i':
				input_filename = optarg;
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 't':
				tau = atoi(optarg);
				break;
			default:
				printf("Option incorrect\n");
				exit(1);
		}
	}
	FILE *fp = fopen(input_filename, "r");
	double qk = 0, qk_sq = 0, qk_qkplustau = 0;
	double *buf = malloc( (tau+1) * sizeof(double));
	int i = 0;
	while(fscanf(fp, "%lf\n", &buf[i % (tau+1)]) == 1) {
		qk += buf[i % (tau+1)];
		qk_sq += buf[i % (tau+1)] * buf[i % (tau+1)];
		if(i >= tau)
			qk_qkplustau += buf[i % (tau+1)] * buf[ (i-tau) % (tau+1)];
		i++;
	}
	fclose(fp);
	free(buf);
	qk /= i;
	qk_sq /= i;
	qk_qkplustau /= i - tau;
	
	double autocorrelation = (qk_qkplustau - qk * qk) / (qk_sq - qk * qk);
	fp = fopen(output_filename, "a");
	fprintf(fp, "%d %lf\n", tau, autocorrelation);
	fclose(fp);
	return 0;

}
