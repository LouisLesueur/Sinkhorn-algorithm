#include "simplex.h"
#include "barycenter.h"
#include "utils.h"
#include <getopt.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
using namespace std;


//Default parameters
int n_iter = 10, steps = 1;
float lambda = 0.5, lambdai=0, lambdao=1, fact = 1000;
string in="../images/circle.png", out="../images/square.png";

void PrintHelp()
{
	cout <<
		"--in [path] <i>:              Set the starting image\n"
		"--out [path] <o>:             Set the final image\n"
		"--n_iter [number] <n>:          Set number of iterations\n"
		"--fact [number] <f>            Set the precision\n"
		"--lambda [number] <l>:          Set the value of lambda\n"
		"--lambdai [number] <I>:        Set the start value of lambda for a serie\n"
		"--lambdao [number] <O>:        Set the stop value of lambda for a serie\n"
		"--steps [number] <s>:           Set the number of steps\n"
		"--help:                Show help\n";
	exit(1);
}

void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "i:o:n:f:l:I:O:s:h";
	const option long_opts[] = {
		{"in", required_argument, nullptr, 'i'},
		{"out", required_argument, nullptr, 'o'},
		{"n_iter", required_argument, nullptr, 'n'},
		{"fact", required_argument, nullptr, 'f'},
		{"lambda", required_argument, nullptr, 'l'},
		{"lambdai", required_argument, nullptr, 'I'},
		{"lambdao", required_argument, nullptr, 'O'},
		{"steps", required_argument, nullptr, 's'},
		{"help", no_argument, nullptr, 'h'},
		{nullptr, no_argument, nullptr, 0}
	};
	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
		if (-1 == opt){
			cout<<"Using default arguments !"<<endl;
			break;
		}
		switch (opt)
		{
			case 'i':
				in = "../images/" + string(optarg);
				cout<<"Starting from: "<<in<<endl;
				break;
			case 'o':
				out = "../images/" + string(optarg);
				cout<<"Going to: "<<out<<endl;
				break;
			case 'n':
				n_iter = stoi(optarg);
				cout << "n_iter set to: " << n_iter << endl;
				break;
			case 'f':
				fact = stof(optarg);
				cout << "fact set to: "<<fact<<endl;
				break;
			case 'l':
				lambda = stof(optarg);
				cout << "lambda set to: " << lambda << endl;
				break;
			case 'I':
				lambdai = stof(optarg);
				cout << "lambdai set to: " << lambdai << endl;
				break;
			case 'O':
				lambdao = stof(optarg);
				cout << "lambdao set to: " << lambdao << endl;
				break;
			case 's':
				steps = stoi(optarg);
				cout <<steps<< " steps" << endl;
				break;

			case 'h': // -h or --help
			case '?': // Unrecognized option
			default:
				PrintHelp();
				break;
		}
	}
}



int main(int argc, char **argv)
{
	ProcessArgs(argc, argv);
	simplex IN(in);
	simplex OUT(out);
	matrix<float> K;

	float eps=1/float(fact*IN.length());

	cout<<"Generating K..."<<endl;
	K = gen_K(IN.length(), eps);

	if(steps==1){
		    cout<<"Building bary.png for n_iter ="<<n_iter<<" eps ="<<eps<<" and lambda ="<<lambda<<endl;
				simplex barycenter = bar(K, IN, OUT, lambda, eps, n_iter, "../out/bary.png");
				barycenter.export_to_img();
				return 0;
	}
	else
	{
		lambda = lambdai;

		for(int i=0; i<steps; i++){
			lambda = (1-(float(i)/float(steps-1)))*lambdai + (float(i)/float(steps-1))*lambdao;
			string name = "../out/bary";
			name += to_string(i);
			name += ".png";
			cout<<"building "<<name<<" for n_iter ="<<n_iter<<" eps ="<<eps<<" and lambda ="<<lambda<<endl;
			simplex barycenter = bar(K, IN, OUT, lambda, eps, n_iter, name);
			cout<<endl;
			barycenter.export_to_img();
		}
	}
	return 0;
}
