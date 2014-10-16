

#include <cmath>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>

using namespace std;


// compute the L2 Norm of the vector to find out its appropriate magnitude
double l2_norm(const vector<double>& v1)
{

	double running_total =0.0;
	for_each(v1.begin(),v1.end(),[&](double val){running_total += (val * val);});
	return sqrt(running_total); // returns the L2 Norm of the vector v1, by sqrt(sum of squares)
}
 
// performs operation Op on both vectors, and stores the results in a temporary vector
template<typename Op>
vector<double> vec_operate(const vector<double>& v1, const vector<double>& v2,Op todo)
{
	vector<double> resultant; // stores the results
	for(int i =0; i < v1.size(); ++i)
	{
		resultant.push_back(todo(v1[i],v2[i])); // stores the results of todo() called on the i'th element on each vector
	}

	return resultant;
} 

// count the in-degree of a particular node, this is the number of elements in column node > 0
int get_in_degree(int node, const vector<vector<double>>& M)
{
	int in = 0;
	for(int i =0; i < M.size(); ++i)
	{
		in += (M[node][i] != 0 ? 1 : 0);
	}
	return in;
}
// count the out-degree of a particular node, this is the number of elements in row node > 0
int get_out_degree(int node,const vector<vector<double>>& M)
{
	int out = 0;
	for(int i =0; i < M.size(); ++i)
	{
		out += (M[i][node] != 0 ? 1 : 0);
	}
	return  out;
}

vector<double> pagerank_converge(const vector<vector<double>>& M,int N, double epsilon,double init, double sum_constraint, double beta)
{
	double diff = 1 << 10; // make absurdly high for now for the while loop

	vector<double> R(N); // page-rank vector


	fill(R.begin(),R.end(), init);

	
	while(diff > epsilon) // run until it converges
	{
		vector<double> r_new(N); // create a new temporary vector to hold the new page-rank values
		for(int i =0; i < N; ++i)
		{
			//cout << i << " = ";
			if(!get_in_degree(i,M)) // if node i doesn't have any incoming edges, no need to calculate its rank (its 0)
			{
			//	cout << "0\n\n";
				r_new[i] = 0.0;
				continue;
			}
			
			for(int j =0 ; j < N; ++j)
			{
				if(M[i][j] != 0) // check if there's a edge from j --> i
				{
					double a_val =  R[j] / get_out_degree(j,M); // compute the partial sum
					//cout << beta << " * " << j << " ( " << R[j] << " ) / " << get_out_degree(j,M) << "\n";
					r_new[i] += beta * a_val; // factor in beta
					//cout << " + ";
				}
			}
			//cout << "\n\n";
		}
		double s = accumulate(r_new.begin(),r_new.end(),0.0); // accumulate the new values

		vector<double> to_add(N); // we need to add a 3x1 vector, whose elements are 1-S/N
		fill(to_add.begin(),to_add.end(),  (sum_constraint-s)/N);

		// return a vector whose elements are the result of r_new[i] + to_add[i] for all i
		// this is the same as vector addition
		r_new = vec_operate(r_new,to_add, [](double i, double j) -> double {return i+j;}); 

		// return a vector whose elements are the result of r_new[i]-R[i]
		// this is used to compute the difference. Once the L2-Norm of the difference is smaller than 
		// epsilon, we've converged
		diff = l2_norm(vec_operate(r_new,R,[](double i, double j) -> double {return i-j;}));
		R = r_new;


	}
	return R;
	
}

vector<double> pagerank_iterate(const vector<vector<double>>& M, int N,double init,double sum_constraint, double beta, int iter)
{
	double diff = 1 << 10; // make absurdly high for now for the while loop

	vector<double> R(N); // page-rank vector
	fill(R.begin(),R.end(),init);

	for(int iter_num =0 ; iter_num < iter; ++iter_num)
	{
		vector<double> r_new(N); // create a new temporary vector to hold the new page-rank values
		for(int i =0; i < N; ++i)
		{
			if(!get_in_degree(i,M)) // if node i doesn't have any incoming edges, no need to calculate its rank (its 0)
			{
				r_new[i] = 0.0;
				continue;
			}
			
			for(int j =0 ; j < N; ++j)
			{
				if(M[i][j] != 0) // check if there's a edge from j --> i
				{
					double a_val =  R[j] / get_out_degree(j,M); // compute the partial sum
					r_new[i] += beta * a_val; // factor in beta
				}
			}
		}
		double s = accumulate(r_new.begin(),r_new.end(),0.0); // accumulate the new values

		vector<double> to_add(N); // we need to add a Nx1 vector, whose elements are 1-S/N
		fill(to_add.begin(),to_add.end(),  (sum_constraint-s)/N);

		// return a vector whose elements are the result of r_new[i] + to_add[i] for all i
		// this is the same as vector addition
		R = vec_operate(r_new,to_add, [](double i, double j) -> double {return i+j;}); 
	}

	return R;
}
int main()
{
	cout << "Number of nodes in the graph\n\n";
	int N; // number of nodes in the graph
	cin >> N;
	cout << "Please enter " << N << " lines of info with the format SOURCE_NODE NUMBER_OF_OUTGOING_EDGES DESTINATION_1 DESTINATION_2 etc..\n\n";
	std::vector<std::vector<double>> M(N); 
	
	// fill in each vector with N slots
	for_each(M.begin(),M.end(),[&](std::vector<double>& vec) { vec.resize(N);});



	for(int i =0 ; i < N; ++i)
	{
		// read in the graph in the following format:
		// SRC NUM_EDGES DST1 DST2 DST3..... DSTN
		int cur_src, num_edges, cur_dst;
		cin >> cur_src >> num_edges;
		for(int j =0 ; j < num_edges; ++j)
		{
			cin >> cur_dst;
			M[cur_dst][cur_src] = 1.0 / num_edges;
		}
	}

	int iterations = 0;
	double epsilon = 0.0;
	double initial = 0.0;
	double sum_to = 0.0; 
	double beta = 0.0;
	cout << "Enter your preferred epsilon value for the convergence method\n\n";
	cin >> epsilon;
	cout << "Enter the amount of iterations to perform for the iterated method\n\n";
	cin >> iterations;
	cout << "Enter the initial values for the page-rank vector (0 for default)\n\n";
	cin >> initial;
	cout << "Enter the summation constraint\n\n";
	cin >> sum_to;
	cout << "Enter the beta-teleport value\n\n";
	cin >> beta;

	initial = initial == 0.0 ? 1.0 / N : initial;

	auto pagerank_convergence  = pagerank_converge(M,N,epsilon,initial,sum_to,beta);//pagerank_converge(M,N,epsilon,beta);
	auto pagerank_iteration =    pagerank_iterate(M,N,initial,sum_to,beta,iterations); 
	
	cout << "page rank until " << epsilon << " convergence\n";
	for_each(pagerank_convergence.begin(),pagerank_convergence.end(),[](double x){cout << x << "\n";});
	cout << endl << "------------------------" <<  endl;
	cout << "Page rank after " << iterations << " iterations\n";
	for_each(pagerank_iteration.begin(),pagerank_iteration.end(),[](double x){cout << x << "\n";});

}