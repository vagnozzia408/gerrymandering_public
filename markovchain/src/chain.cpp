#include <iostream>
#include <iomanip>
//#include <math.h>
#include <cmath>
#include <string.h>
#include <stdio.h>
#include <random>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "chaincmdline.h"
#include <queue>
#include <unordered_map>

//#define NDEBUG

#define MAXDEGREE 50 //max number of precincts incident with any precinct.
using namespace std;



enum Type_S {
	rootep=0,
	GCserial=1,
	starsplit=2,
	GCparallel=3
};

//GLOBAL VARIABLES
int g_MAXEDGES;      //theoretical maximum districting degree
int g_NUMDISTRICTS;  //number of districts in the state
char g_debuglevel;   //debug output level
Type_S StatTest;     //which statistical test do we use?
int g_seedoffset;    //random seed offset
////////////////////////



struct edge{
	int u;
	int j;
	edge(int a, int b){
		u=a;
		j=b;
	}
	edge(){}
	operator int () {   //assigns an unique integer to an edge
		return MAXDEGREE*this->u+j;
	}
};


void Tokenize( const char* szString, vector<std::string>&
vecstrTokens, const char* szSeparators, bool fNoEmpties){
	const char* pc;
	string strCur;
	bool fPush;

	if( !( pc = szString ) )
	return;

	fPush = false;
	while( true ) {
		strCur.clear( );
		if( fNoEmpties )
		for( ; *pc && strchr( szSeparators, *pc ); ++pc );
		if( !*pc ) {
			if( !fNoEmpties && fPush )
			vecstrTokens.push_back( strCur );
			return; }
		for( ; *pc && !strchr( szSeparators, *pc ); ++pc )
		strCur += *pc;
		if( (fPush = !!*pc) )
		pc++;
		vecstrTokens.push_back( strCur ); 
	} 
}


class precinct{
	bool initialized;
public:
	bool selfinitialized;
	int degree;
	int * neighbors; //list of INDICES (in pr[]) of neighbors of the precinct
	int * self; //my place in each of my neighbor's neighbor lists
	double * shared_perimeters; //list of shared perimeters with the neighbors
	double area;
	int population;
	int voteA;
	int voteB;
	int original_district;
	string county;
	bool frozen;
	bool giant; //Am I in the intial component of the original district.
	friend void computeselves();

	precinct(){
		initialized=false;
		selfinitialized=false;
		giant=false;
		frozen=false;
	}

	precinct(int d){
		degree=d;
		neighbors=new int[degree];
		self=new int[degree];
		shared_perimeters= new double[degree];
		initialized=selfinitialized=true;
		frozen=false;
		giant=false;
	}

	precinct(vector<string> dataline, bool flip, bool use_counties){
		initialized=true;
		selfinitialized=false;
		vector<string> neighbors_string;
		vector<string> shared_perimeters_string;
		Tokenize(dataline[1].c_str(),neighbors_string,",",true);
		Tokenize(dataline[2].c_str(),shared_perimeters_string,",",true);
		degree=neighbors_string.size();
		if (shared_perimeters_string.size()!=degree){
			cerr << "ERROR: perimeter list differs in length from neighbor list?" << endl;
			cerr << "neighbor list was: "<< dataline[1] << endl; 
			exit(-1);
		}
		neighbors=new int[degree];
		shared_perimeters= new double[degree];
		
		for (int i=0; i<degree; i++){
			neighbors[i]=atoi(neighbors_string[i].c_str()); //subtract one because of indexing from 0
			shared_perimeters[i]=atof(shared_perimeters_string[i].c_str());
		}
		area=atof(dataline[4].c_str());
		population=atoi(dataline[5].c_str());
		if (flip){
			voteA=atoi(dataline[7].c_str());
			voteB=atoi(dataline[6].c_str());
		}
		else{
			voteA=atoi(dataline[6].c_str());
			voteB=atoi(dataline[7].c_str());
		}
		
		original_district=atoi(dataline[8].c_str())-1; //subtract one because of indexing from 0
		if (use_counties)
		county=dataline[9];
		frozen=false;
		giant=false;
		if (original_district>=g_NUMDISTRICTS){
			cerr << "ERROR: district number out of range"<<endl;
			exit(-1);
		}
	}

	precinct(const precinct& source){
		initialized=source.initialized;
		selfinitialized=source.selfinitialized;
		degree=source.degree;
		area=source.area;
		population=source.population;
		voteA=source.voteA;
		voteB=source.voteB;
		original_district=source.original_district;
		county=source.county;
		frozen=source.frozen;
		giant=source.giant;
		if (initialized){
			neighbors=new int[degree];
			shared_perimeters=new double[degree];
			for (int i=0; i<degree; i++){
				neighbors[i]=source.neighbors[i];
				shared_perimeters[i]=source.shared_perimeters[i];
			}
		}
		if (selfinitialized){
			self=new int[degree];
			for (int i=0; i<degree; i++)
			self[i]=source.self[i];
		}
	}
	precinct& operator= (const precinct& source){
		if (&source!=this){
			initialized=source.initialized;
			selfinitialized=source.selfinitialized;
			degree=source.degree;
			area=source.area;
			population=source.population;
			voteA=source.voteA;
			voteB=source.voteB;
			original_district=source.original_district;
			county=source.county;
			frozen=source.frozen;
			giant=source.giant;
			if (initialized){
				neighbors=new int[degree];
				shared_perimeters=new double[degree];
				for (int i=0; i<degree; i++){
					neighbors[i]=source.neighbors[i];
					shared_perimeters[i]=source.shared_perimeters[i];
				}
			}
			if (selfinitialized){
				self=new int[degree];
				for (int i=0; i<degree; i++)
				self[i]=source.self[i];
			}
		}
		return *this;
	}
	~precinct(){
		if (initialized){
			delete[] neighbors;
			delete[] shared_perimeters;
		}
		if (selfinitialized)
		delete[] self;
	}
};


void computeselves(precinct * pr, int N){   //find edges back to me 
	//from my neighboring precints.
	//where multiple edges are present between a pair,
	//we don't need or ensure geometric correspondence.  Only 1-1 correspondence.

	int * boundary;//indices of precincts bordering boundary, in clockwise order
	int * bself; //boundary back pointers
	boundary=new int[N];
	bself=new int[N];

	for (int i=0; i<N; i++){
		boundary[i]=-1;  
		bself[i]=-1;
	}

	int boundarydegree=0;

	for (int i=0; i<N; i++){
		pr[i].selfinitialized=true;
		pr[i].self=new int[pr[i].degree];
		for (int j=0; j<pr[i].degree; j++){
			if (pr[i].neighbors[j]>-1){
				int setcount=0;
				for (int l=0; l<pr[pr[i].neighbors[j]].degree; l++){
					//if (pr[pr[i].neighbors[j]].neighbors[l]==i && pr[pr[i].neighbors[j]].shared_perimeters[l]==pr[i].shared_perimeters[j]){
					// AMV EDIT: Adjusted statement to include a tolerance for shared perimeter lengths vs. strict equality. This accounts for rounding error during neighbor calculations.
					if (pr[pr[i].neighbors[j]].neighbors[l]==i && abs(pr[pr[i].neighbors[j]].shared_perimeters[l]-pr[i].shared_perimeters[j])<1e-10){
						setcount++;
						pr[i].self[j]=l;
					}
				}
				assert(setcount<2);
				if (setcount!=1){
					cerr << "setcount is " << setcount<<endl;
					cerr << "i is "<< i<<endl;
					cerr << "j is "<< j<<endl;
					cerr << "pr[i].neighbors[j] is " << pr[i].neighbors[j] << endl;
				}
				assert(setcount==1);
			}
			else{
				// QUESTION!!! Why did I add an if statement here?
				if (boundary[-pr[i].neighbors[j]-1]==-1){
					boundarydegree++;
					boundary[-pr[i].neighbors[j]-1]=i;  
					bself[-pr[i].neighbors[j]-1]=j;

					pr[i].self[j]= (-pr[i].neighbors[j]-1); //outside edges are labeled by order
					pr[i].neighbors[j]=N; //reset pointer to boundary precinct
				}
				else{
					cerr << "ERROR 276: Precinct i=" << i << endl;
				}
			}
		}
	}
	
	assert(boundary[boundarydegree]<0);
	assert(boundary[boundarydegree-1]>=0);

	assert(bself[boundarydegree]<0);
	assert(bself[boundarydegree-1]>=0);


	pr[N]=precinct(boundarydegree);
	for (int j=0; j<pr[N].degree; j++){
		pr[N].neighbors[j]=boundary[j];
		pr[N].self[j]=bself[j];
		pr[N].shared_perimeters[j]=pr[pr[N].neighbors[j]].shared_perimeters[ pr[N].self[j] ];
	}
	pr[N].original_district=-1; //dummy district
	
	delete [] boundary;
	delete [] bself;
}

//CODE to check for counties which belong to only one district and set
//their precincts to frozen:
void freezedistrictsbycounty(precinct * pr, int N){   
	unordered_map<string,int> countyfreeze_map; //unique dists for frozen counties
	for (int i=0; i<N; i++){
		if (countyfreeze_map.count(pr[i].county)==0)
		countyfreeze_map[pr[i].county]=pr[i].original_district;
		else if (countyfreeze_map[pr[i].county]!=pr[i].original_district) 
		countyfreeze_map[pr[i].county]=-1; //not preserved
	}
	
	for (auto keyvaluepair : countyfreeze_map){
		//    cout << "first is " << keyvaluepair.first<< " and second is "<<keyvaluepair.second<<endl;
		if (keyvaluepair.second!=-1)
		cout << keyvaluepair.first << " County is preserved by initial districting"<<endl;
	}
	
	for (int i=0; i<N; i++){
		if (countyfreeze_map[pr[i].county]!=-1){
			assert(countyfreeze_map[pr[i].county]==pr[i].original_district);
			pr[i].frozen=true;
		}
		else
		assert(pr[i].frozen==false);
	}
}



template<class TYPE>  //type must be uniquely castable to (int) 
class rdpile{ //fast deletion of specific elements and fast sampling of random elements
private:
	bool initialized;
	TYPE * array;
	int max; //max number of TYPES to be stored
	int maxindex; //max size of (int) TYPE to be seen
	int * index; //index of locations
public:
	int count; //number of elements

	TYPE access(int i){
		assert (i>=0 && i<count);
		return array[i];
	}
	void insert(TYPE e){
		index[(int) e]=count;
		array[count++]=e;
	}
	void removeindex(int i){
		assert (i>=0 && i<count);
		array[i]=array[--count];
		index[(int) array[i]]=i;
	}
	void remove(TYPE e){
		int i=index[(int) e];
		assert (i>=0 && i<count);
		array[i]=array[--count];
		index[(int) array[i]]=i;
	}
	rdpile(){
		initialized=false;
	}
	rdpile(int m, int M){
		initialized=true;
		max=m;
		maxindex=M;
		array=new TYPE[max];
		index=new int[maxindex];
		count=0;
	}
	rdpile(const rdpile& source){
		if (source.initialized==true){
			initialized=true;
			max=source.max;
			maxindex=source.maxindex;
			count=source.count;
			array=new TYPE[source.max];
			index=new int[source.maxindex];
			for (int i=0; i<count; i++){
				array[i]=source.array[i];
				index[(int) array[i]]=source.index[(int) array[i]];
			}
		}
		else
		initialized=false;
	}
	rdpile& operator= (const rdpile& source){
		if (&source!=this){
			if (source.initialized==true){
				initialized=true;
				max=source.max;
				maxindex=source.maxindex;
				count=source.count;
				array=new TYPE[source.max];
				index=new int[source.maxindex];
				for (int i=0; i<count; i++){
					array[i]=source.array[i];
					index[(int) array[i]]=source.index[(int) array[i]];
				}
			}
			else
			initialized=false;
		}
		return *this;
	}
	~rdpile(){
		if (initialized){
			delete[] array;
			delete[] index;
		}
	}
};



bool validpop(int pop, double avgpop, double popthresh){
	return pop<avgpop*(1+popthresh) && pop>avgpop*(1-popthresh);
}


double dcompact(double Area, double Perim){  //inverse of Polsby-Popper: bigger is worse
	return (pow(Perim,2)/(Area*4*M_PI));
}

double compactsum(double * Area, double * Perim, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=dcompact(Area[i],Perim[i]);
	return sum;
}

double compactsum(int * Area, double * Perim, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=dcompact(Area[i],Perim[i]);
	return sum;
}

double compactL2(double * Area, double * Perim, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=pow(dcompact(Area[i],Perim[i]),2);
	return pow(sum,.5);
}

double compactL2(int * Area, double * Perim, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=pow(dcompact(Area[i],Perim[i]),2);
	return pow(sum,.5);
}

double arraysum(double * array, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=array[i];
	return sum;
}


double arrayL2(double * array, int arraylength){
	double sum=0;
	for (int i=0; i<arraylength; i++)
	sum+=pow(array[i],2);
	return pow(sum,.5);
}


int reps(double * Ashare){
	int reps=0;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		if (Ashare[k]>.5)
		reps++;
	}
	return reps;
}

double variance(double * Ashare)
{
	double meansq=0;
	double mean=0;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		meansq+=pow(Ashare[k],2);
		mean+=Ashare[k];
	}
	meansq=meansq/g_NUMDISTRICTS;
	mean=mean/g_NUMDISTRICTS;
	return meansq-pow(mean,2);
}

double efficiency_gap(int * DvotesA, int * DvotesB)
{
	int total=0;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		total+=DvotesA[k];
		total+=DvotesB[k];
	}
	int wastedA=0;
	int wastedB=0;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		if (DvotesA[k]>DvotesB[k]){
			wastedA+=(DvotesA[k]-DvotesB[k])/2;
			wastedB+=DvotesB[k];
		}
		else{
			wastedB+=(DvotesB[k]-DvotesA[k])/2;
			wastedA+=DvotesA[k];
		}
	}
	return ((double) wastedA-wastedB)/total;
}

double median_mean(double * Ashare)
{
	double * AshareCpy;
	AshareCpy = new double[g_NUMDISTRICTS];
	copy(Ashare, Ashare+g_NUMDISTRICTS,AshareCpy);
	double mean,median,sum;
	int middle=(g_NUMDISTRICTS+1)/2; //top median
	nth_element(AshareCpy, AshareCpy+middle, AshareCpy+g_NUMDISTRICTS); // Sorts elements until the nth element is in the nth position
	nth_element(AshareCpy, AshareCpy+middle-1, AshareCpy+middle);		// In other words, it doesn't sort the whole thing
	if (!(g_NUMDISTRICTS%2))//for median with even g_NUMDISTRICTS
		median=(AshareCpy[middle]+AshareCpy[middle-1])/2;
	else  //median for odd g_NUMDISTRICTS
		median=AshareCpy[middle-1];
	sum=0;
	for (int k=0; k<g_NUMDISTRICTS; k++)
		sum+=AshareCpy[k];
	mean=sum/g_NUMDISTRICTS;
	delete [] AshareCpy;
	return mean-median;
}

////////////////////////////////////////////////////////////////////////////////////

// AMV EDIT #1: Added a function to compute the bias for a particular districting.

double partisan_bias(double * Ashare)
{
	// Calculate the average vote share for party A under the current districting.
	double sum_Ashare = 0;
	// Calculate the estimated seat share for party A under the current districting.
	int sum_seats = 0;
	for (int i=0; i<g_NUMDISTRICTS; i++){
		sum_Ashare += Ashare[i];
		if (Ashare[i] > 0.5){
			sum_seats += 1;
		}
	}
	double avg_Ashare = sum_Ashare / g_NUMDISTRICTS;	
	double seat_share = double(sum_seats) / g_NUMDISTRICTS; // Cast sum_seats to double so division works.
	
	// What if the parties switched? If Party B won the same statewide proportion of seats, would the estimated seat count change?
	// Define the uniform partisan swing applied to each district vote proportion for Party A.
	double swing = 1-2*avg_Ashare;
	
	double * Ashare_new;
	int new_seats = 0;
	Ashare_new = new double[g_NUMDISTRICTS];
	
	for (int i=0; i<g_NUMDISTRICTS; i++){
		Ashare_new[i] = Ashare[i]+swing;
		if (Ashare_new[i] > 0.5){
			new_seats += 1;
		}
	}
	double new_seat_share = double(new_seats) / g_NUMDISTRICTS;
		
	// Estimate the bias measure.
	double bias	= (seat_share-(1-new_seat_share))/2;
	delete [] Ashare_new;	
	return bias;
}

////////////////////////////////////////////////////////////////////////////////////////////////

// AMV EDIT 10-26-19: Compute Geometric Partisan Bias BG under Uniform Partisan Swing Assumption

double BG_uniform(double * Ashare)
{
	// Ashare: Democratic Vote Shares per District
	// Calculate the average vote share V for party A under the current districting.
	// Calculate the estimated seat share S(V) for party A under the current districting.
	double sum_Ashare = 0;
	int sum_seats = 0;
	for (int i=0; i<g_NUMDISTRICTS; i++){
		sum_Ashare += Ashare[i];
		if (Ashare[i] > 0.5){
			sum_seats += 1;
		}
	}
	double V = sum_Ashare / g_NUMDISTRICTS; 
	double SV = double(sum_seats) / g_NUMDISTRICTS; // Cast sum_seats to double so division works.
	
	// Create arrays to store seats-votes curve points. 
	static double * V_Points = new double[2*g_NUMDISTRICTS];
	static double * SV_Points = new double[2*g_NUMDISTRICTS];
	
	// Store original point.
	V_Points[0] = V;
	SV_Points[0] = SV;
	
	// Calculate seats-votes curve points.
	static double * Ashare_new = new double[g_NUMDISTRICTS];	// Each iteration overwrites the old stuff
	
	for (int i=0; i<g_NUMDISTRICTS; i++){				// for each district...
		// First calculate the partisan swing needed to flip this seat.
		double delta = 0;
		if (Ashare[i] < 0.5){					// TO DO: THE EQUALITY SHOULD PROBABLY BE HERE. OK.
			delta = 0.50001 - Ashare[i];
		}
		else { // if (Ashare[i] >= 0.5) (equality seldom occurs)
			delta = 0.49999 - Ashare[i];
		}
		
		// Calculate new district shares.
		// Calculate new V (avg_Ashare_new) and add to V_Points.
		// Calculate new SV (seat_share_new) and add to SV_Points.
		double sum_Ashare_new = 0;
		int sum_seats_new = 0;
		for (int j=0; j<g_NUMDISTRICTS; j++){
			Ashare_new[j] = Ashare[j]+delta;
			sum_Ashare_new += Ashare_new[j];
			if (Ashare_new[j] > 0.5){
				sum_seats_new += 1;
			}
		}
		double newV = sum_Ashare_new / g_NUMDISTRICTS;
		V_Points[i+1] = newV;
		double newSV = double(sum_seats_new) / g_NUMDISTRICTS;
		SV_Points[i+1] = newSV;
	}
	
	// Sort the points stored so far. This is okay since the seats-votes curve is strictly increasing.
	sort(V_Points, V_Points + g_NUMDISTRICTS + 1);
	sort(SV_Points, SV_Points + g_NUMDISTRICTS + 1);
	
	// Create arrays to store the inverted seats-votes curve points.
	static double * IV_Points = new double[2*g_NUMDISTRICTS];
	static double * ISV_Points = new double[2*g_NUMDISTRICTS];
	for (int i=0; i<g_NUMDISTRICTS+1; i++){
		IV_Points[i] = 1 - V_Points[g_NUMDISTRICTS-i];
		ISV_Points[i] = 1 - SV_Points[g_NUMDISTRICTS-i];
	}
	
	// Find the intersection points.
	int k=0;
	for (int i=0; i<g_NUMDISTRICTS; i++){
		if ( ((V_Points[i] < IV_Points[i]) && (V_Points[i+1] > IV_Points[i+1])) || ((V_Points[i] > IV_Points[i]) && (V_Points[i+1] < IV_Points[i+1])) ){
			// Seats-Votes Line Segment
			double m1 = (SV_Points[i+1] - SV_Points[i])/(V_Points[i+1]-V_Points[i]);
			double b1 = SV_Points[i] - m1*V_Points[i];
			
			// Inverted Seats-Votes Line Segment
			double m2 = (ISV_Points[i+1] - ISV_Points[i])/(IV_Points[i+1]-IV_Points[i]);
			double b2 = ISV_Points[i] - m2*IV_Points[i];
			
			// Intersection Point
			double x = (b2-b1)/(m1-m2);
			double y = m1 * x + b1;
			
			// Add the intersection point. 
			V_Points[g_NUMDISTRICTS+k+1] = x;
			IV_Points[g_NUMDISTRICTS+k+1] = x;
			SV_Points[g_NUMDISTRICTS+k+1] = y;
			ISV_Points[g_NUMDISTRICTS+k+1] = y;
			
			k++;
		}
	}
	
	// Sort all the points.
	sort(V_Points, V_Points + g_NUMDISTRICTS + k + 1);
	sort(SV_Points, SV_Points + g_NUMDISTRICTS + k + 1);
	sort(IV_Points, IV_Points + g_NUMDISTRICTS + k + 1);
	sort(ISV_Points, ISV_Points + g_NUMDISTRICTS + k + 1);
	
	// Find area between the two seats-votes curves.	
	double uniform_geom_bias = 0;
	double xmax = max(V_Points[g_NUMDISTRICTS+k],IV_Points[g_NUMDISTRICTS+k]);
	
	for(int i=0; i<g_NUMDISTRICTS+k; i++){
		// Area under Seats-Votes Curve using trapezoids.
		double b1 = xmax - V_Points[i];
		double b2 = xmax - V_Points[i+1];
		double h = SV_Points[i+1] - SV_Points[i];
		double area1 = 0.5 * (b1+b2) * h;
		
		// Area under Inverted Seats-Votes Curve.
		double ib1 = xmax - IV_Points[i];
		double ib2 = xmax - IV_Points[i+1];
		double ih = ISV_Points[i+1] - ISV_Points[i];
		double area2 = 0.5 * (ib1+ib2) * ih;
		
		uniform_geom_bias += abs(area2-area1);
	}
	
	return uniform_geom_bias;
}
////////////////////////////////////////////////////////////////////////////////////////////////

// AMV EDIT 11-1-19: Compute Geometric Partisan Bias BG under Modified Partisan Swing Assumption

double BG_modified(double * Ashare)
{
	// Ashare: Democratic Vote Shares per District	
	// Calculate the average vote share V for party A under the current districting.
	double sum_Ashare = 0;						// Incorporated randomness to bump any v=.5 up or down slightly
	for (int i=0; i<g_NUMDISTRICTS; i++){
		sum_Ashare += Ashare[i];
	}
	double V = sum_Ashare / g_NUMDISTRICTS; 
	
	static double * V_Points = new double[2*g_NUMDISTRICTS];
	
	// Add the initial average vote share V to the list of average vote shares.
	V_Points[0]=V;
	
	// Calculate new average district vote shares under modified partisan swing.
	for (int i=0; i<g_NUMDISTRICTS; i++){
		double new_V = 0;
		// If Party B occupies the seat...
		if (Ashare[i] < 0.5){
			// The seat will be lost if the statewide vote falls to new_V.
			new_V = 1 - (1-V)/(2*(1-Ashare[i]));
		}
		// If Party A occupies the seat...
		else { // if (Ashare[i] >= 0.5) (equality seldom occurs)		BUT WHEN IT DOES IT IS BAD. WHY DID YOU DO THIS.
			new_V = V / (2*Ashare[i]);									// 3-2-2020 EDITS: Fixed Ashare calculations to ensure no ties result.
		}
		if(isnan(new_V)){ // 2021-11-29: -nan getting calculated as a V_Point. This should locate the error.
			cerr << "V=" << setprecision(20) << V << endl;
			cerr << "Ashare[" << i << "]=" << setprecision(20)<< Ashare[i]  << endl;
			exit(-1);
		}
		
		// 2020-03-02
		//if(new_V == V){
		//	cerr << "new_V == V? Ashare is " << Ashare[i] << endl;
		//	exit(-1);
		//}
		
		// 2020-02-24		
		//if (isnan(new_V)){
		//	cerr << "new_V is " << new_V << endl;
		//	cerr << "i=" << i << endl;
		//	exit(-1);
		//}
		/////////////
		
		V_Points[i+1] = new_V;
		// 2021-12-05 Check to see if somehow new_V magically became a nan.
		if(isnan(V_Points[i+1])){
			cerr << "V_Points[" << i+1 << "]=" << setprecision(20) << V_Points[i+1] << endl;
			cerr << "V=" << setprecision(20) << V << endl;
			cerr << "Ashare[" << i << "]=" << setprecision(20)<< Ashare[i]  << endl;
			exit(-1);
		}
	}
	
	// Sort the average vote shares.
	sort(V_Points, V_Points + g_NUMDISTRICTS + 1);
	// Calculate corresponding seat shares, which will already be sorted.
	static double * SV_Points = new double[2*g_NUMDISTRICTS];
	for (int i=0; i<g_NUMDISTRICTS+1; i++){
		SV_Points[i] = double(i)/g_NUMDISTRICTS; // Cast number of seats i to double.
	}
	
	// AMV EDITS 2020-03-02 //
	// Make sure we don't have consecutive points. This also results in division by zero.
	for (int i=0; i<g_NUMDISTRICTS; i++){	// for each pair of V_Points
		if (V_Points[i] == V_Points[i+1]){	// if we observe two consecutive V_Points
			if(i==g_NUMDISTRICTS-1){		// if the last two points are consecutive
				double m = (SV_Points[i+1]-SV_Points[i-1])/(V_Points[i+1]-V_Points[i-1]);
				double b = SV_Points[i-1]-m*V_Points[i-1];
				double adj_V = (SV_Points[i]-b)/m;
				//assert(adj_V > V_Points[i-1] && adj_V < V_Points[i+1]);
				// 11-23-2021 edit: if we hit more than 2 consecutive points, just don't count this measure on this map
				if(adj_V == V_Points[i-1] || adj_V == V_Points[i+1])
					return 2;
				V_Points[i] = adj_V;	// Replace the first consecutive point with this adjusted point.
				// 2021-12-05: Check the replaced V_Point
				if(isnan(V_Points[i])){
					cerr << "First Branch Fails (Line 783)" << endl;
					cerr << "V_Points[" << i-1 << "]=" << setprecision(20) << V_Points[i-1] << endl;
					cerr << "V_Points[" << i << "]=" << setprecision(20) << V_Points[i] << endl;
					cerr << "V_Points[" << i+1 << "]=" << setprecision(20) << V_Points[i+1] << endl;
					cerr << "m=" << setprecision(20) << m << endl;
					cerr << "b=" << setprecision(20) << b << endl;
					cerr << "SV_Points[" << i-1 << "]=" << setprecision(20) << SV_Points[i-1] << endl;
					cerr << "SV_Points[" << i << "]=" << setprecision(20) << SV_Points[i] << endl;
					cerr << "SV_Points[" << i+1 << "]=" << setprecision(20) << SV_Points[i+1] << endl;
					exit(-1);
				}
			}
			else{	// if the two consecutive points are NOT the last two
				double m = (SV_Points[i+2]-SV_Points[i])/(V_Points[i+2]-V_Points[i]);
				double b = SV_Points[i]-m*V_Points[i];
				double adj_V = (SV_Points[i+1]-b)/m;
				//assert(adj_V > V_Points[i] && adj_V < V_Points[i+2]);
				// 11-23-2021 edit: if we hit more than 2 consecutive points, just don't count this measure on this map
				if(adj_V == V_Points[i] || adj_V == V_Points[i+2])
					return 2;
				V_Points[i+1] = adj_V;	// Replace the second consecutive point with this adjusted point.
				// 2021-12-05: Check the replaced V_Point
				if(isnan(V_Points[i+1])){
					cerr << "Second Branch Fails (Line 806)" << endl;
					cerr << "V_Points[" << i << "]=" << setprecision(20) << V_Points[i] << endl;
					cerr << "V_Points[" << i+1 << "]=" << setprecision(20) << V_Points[i+1] << endl;
					cerr << "V_Points[" << i+2 << "]=" << setprecision(20) << V_Points[i+2] << endl;
					cerr << "m=" << setprecision(20) << m << endl;
					cerr << "b=" << setprecision(20) << b << endl;
					cerr << "SV_Points[" << i << "]=" << setprecision(20) << SV_Points[i] << endl;
					cerr << "SV_Points[" << i+1 << "]=" << setprecision(20) << SV_Points[i+1] << endl;
					cerr << "SV_Points[" << i+2 << "]=" << setprecision(20) << SV_Points[i+2] << endl;
					exit(-1);
				}
			}
		}
		// otherwise do nothing - all the V_Points are fine.
	}
	//////////////////////////
	
	// Find inverted seats-votes curve.
	static double * IV_Points = new double[2*g_NUMDISTRICTS];
	static double * ISV_Points = new double[2*g_NUMDISTRICTS];
	for (int i=0; i<g_NUMDISTRICTS+1; i++){					// TO DO: Change for BG_uniform as well - saves a sort call
		IV_Points[i] = 1 - V_Points[g_NUMDISTRICTS-i];
		ISV_Points[i] = 1 - SV_Points[g_NUMDISTRICTS-i];
	}
	
	// Find the intersection points.
	int k=0; // Counts number of intersection points
	for (int i=0; i<g_NUMDISTRICTS; i++){
		
		if ( ((V_Points[i] < IV_Points[i]) && (V_Points[i+1] > IV_Points[i+1])) || ((V_Points[i] > IV_Points[i]) && (V_Points[i+1] < IV_Points[i+1])) ){
			// Seats-Votes Line Segment
			double m1 = (SV_Points[i+1] - SV_Points[i])/(V_Points[i+1]-V_Points[i]);  // DIVISION BY ZERO HAPPENS HERE. How do we get two consecutive V_Points???
			double b1 = SV_Points[i] - m1*V_Points[i];
			
			// Inverted Seats-Votes Line Segment
			double m2 = (ISV_Points[i+1] - ISV_Points[i])/(IV_Points[i+1]-IV_Points[i]);
			double b2 = ISV_Points[i] - m2*IV_Points[i];
			
			// Intersection Point
			double x = (b2-b1)/(m1-m2);
			double y = m1 * x + b1;
			
			// 2020-02-24 
			if(isnan(x) || isnan(y)){
				//cerr << "Ashare[" << i << "] = " << Ashare[i] << " and Ashare[" << i+1 << "] = " << Ashare[i+1] << endl;
				cerr << "SV_Points[" << i << "]=" << setprecision(20) << SV_Points[i] << endl;
				cerr << "SV_Points[" << i+1 << "]=" << setprecision(20) << SV_Points[i+1] << endl;
				cerr << "V_Points[" << i << "]=" << setprecision(20) << V_Points[i] << endl;
				cerr << "V_Points[" << i+1 << "]=" << setprecision(20) << V_Points[i+1] << endl;
				cerr << "Denominator is " << setprecision(20) << V_Points[i+1]-V_Points[i] << endl;
				if (V_Points[i]==V_Points[i+1]){
					cerr << "V_Points are equal." << endl;
				}
				//cerr << "Intersection point k=" << k << " for i=" << i << endl;
				//cerr << "b1=" << b1 << ", b2=" << b2 << ", m1=" << m1 << ", m2=" << m2 << endl;
				//cerr << "x=" << x << " and y=" << y << endl;
				exit(-1);
			}
			/////////////
			
			// Add the intersection point. 
			V_Points[g_NUMDISTRICTS+k+1] = x;
			IV_Points[g_NUMDISTRICTS+k+1] = x;
			SV_Points[g_NUMDISTRICTS+k+1] = y;
			ISV_Points[g_NUMDISTRICTS+k+1] = y;
			
			k++;
		}
	}
	
	// 2020-02-24
	for (int i=0; i<g_NUMDISTRICTS+k; i++){
		if (isnan(V_Points[i]) || isnan(V_Points[i+1]) || isnan(SV_Points[i]) || isnan(SV_Points[i+1])){ // temporarily ignoring IV_Points and ISV_Points
			cerr << "BEFORE sort... (for k = " << k << ")" << endl;
			cerr << "V_Points[" << i << "]=" << V_Points[i] << endl;
			cerr << "V_Points[" << i+1 << "]=" << V_Points[i+1] << endl;
			cerr << "SV_Points[" << i << "]=" << SV_Points[i] << endl;
			cerr << "SV_Points[" << i+1 << "]=" << SV_Points[i+1] << endl;
			exit(-1);
		}
	}
	/////////////
	
	// Sort all the points.
	sort(V_Points, V_Points + g_NUMDISTRICTS + k + 1);
	sort(SV_Points, SV_Points + g_NUMDISTRICTS + k + 1);
	sort(IV_Points, IV_Points + g_NUMDISTRICTS + k + 1);
	sort(ISV_Points, ISV_Points + g_NUMDISTRICTS + k + 1);
	
	// 2020-02-24
	for (int i=0; i<g_NUMDISTRICTS+k; i++){
		if (isnan(V_Points[i]) || isnan(V_Points[i+1]) || isnan(SV_Points[i]) || isnan(SV_Points[i+1])){ // temporarily ignoring IV_Points and ISV_Points
			cerr << "AFTER sort..." << endl;
			cerr << "V_Points[" << i << "]=" << V_Points[i] << endl;
			cerr << "V_Points[" << i+1 << "]=" << V_Points[i+1] << endl;
			cerr << "SV_Points[" << i << "]=" << SV_Points[i] << endl;
			cerr << "SV_Points[" << i+1 << "]=" << SV_Points[i+1] << endl;
			exit(-1);
		}
	}
	/////////////
	
	// Find area between the two seats-votes curves.	
	double modified_geom_bias = 0;
	double xmax = max(V_Points[g_NUMDISTRICTS+k],IV_Points[g_NUMDISTRICTS+k]);
	
	for(int i=0; i<g_NUMDISTRICTS+k; i++){
		// Area under Seats-Votes Curve using trapezoids.
		double b1 = xmax - V_Points[i];
		double b2 = xmax - V_Points[i+1];
		double h = SV_Points[i+1] - SV_Points[i];
		double area1 = 0.5 * (b1+b2) * h;
		
		// Area under Inverted Seats-Votes Curve.
		double ib1 = xmax - IV_Points[i];
		double ib2 = xmax - IV_Points[i+1];
		double ih = ISV_Points[i+1] - ISV_Points[i];
		double area2 = 0.5 * (ib1+ib2) * ih;
		
		// 2020-02-24
		if (isnan(b2) || isnan(h) || isnan(ib2) || isnan(ih)){
			cerr << "b1=" << b1 << ", b2=" << b2 << ", h=" << h << ", ib1=" << ib1 << ", ib2=" << ib2 << ", ih=" << ih << endl;
			cerr << "V_Points[" << i << "]=" << V_Points[i] << endl;
			cerr << "V_Points[" << i+1 << "]=" << V_Points[i+1] << endl;
			cerr << "SV_Points[" << i+1 << "]=" << SV_Points[i+1] << endl;
			cerr << "SV_Points[" << i << "]=" << SV_Points[i] << endl;
		exit(-1);
		}
		/////////////
		
		modified_geom_bias += abs(area2-area1);
		
	}
	
	// 2020-02-24
	if (isnan(modified_geom_bias)){
		cerr << "BG_modified is " << modified_geom_bias << endl;
		exit(-1);
	}
	return modified_geom_bias;
}

////////////////////////////////////////////////////////////////////////////////////

double seat_slide(double * Ashare)
{
	int seats=0;
	double slidedown=.9;
	double slideup=.9;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		if (Ashare[k]>.5){
			seats++;
			slidedown=min(slidedown,Ashare[k]-.5);
		}
		else {
			slideup=min(slideup,.5-Ashare[k]);
		}
	}
	return seats+.5+slidedown-slideup;
}

void tosvg(char* svgfilename, char* inputsvgfilename, precinct* pr, int* currentdistrict, int firstline, int N)
{
	ofstream districtingsvg(svgfilename);
	ifstream inputsvg(inputsvgfilename);
	if (!districtingsvg.good() || !inputsvg.good()){
		cerr << "ERROR with one of the svg files"<<endl;
		exit(-1);
	}
	string correct_header ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");


	mt19937_64 gen(314159265358979+g_seedoffset); //reinitialize so colors are consistent on each run
	std::uniform_int_distribution<int> randpercent(0,100);
	string * color;
	color=new string[g_NUMDISTRICTS];
	for (int i=0; i<g_NUMDISTRICTS; i++){
		int red=randpercent(gen);
		int green=randpercent(gen);
		int blue=randpercent(gen);
		char tmpcolor[25];
		sprintf(tmpcolor,"rgb(%d%%,%d%%,%d%%);",red,blue,green);
		color[i]=string(tmpcolor);
	}
	string line;
	getline (inputsvg,line); //file header
	districtingsvg << line << endl;
	if (correct_header.compare(line)){
		cerr << "ERROR: incorrect file header at pos " << correct_header.compare(line) << endl;
		exit(-1);
	}

	for (int i=0; i<firstline-2; i++){
		getline (inputsvg,line);
		districtingsvg << line << endl;
	}

	for (int i=0; i<N; i++){
		assert(currentdistrict[i]>=0);
		assert(currentdistrict[i]<g_NUMDISTRICTS);
		getline (inputsvg,line); //color line ignored
		districtingsvg << color[currentdistrict[i]]<<endl;
		getline (inputsvg,line); 
		districtingsvg << line << endl;    
	}
	while (getline (inputsvg,line))
	districtingsvg << line;    
	districtingsvg.close();
	inputsvg.close();
}


void tofile(char* filename, precinct* pr, int* currentdistrict, int N, bool use_counties)
{
	ofstream output(filename);
	if (!output.good()){
		cerr << "ERROR with precinct output file"<<endl;
		exit(-1);
	}
	if (use_counties)
	output << "precinctlistv02" <<endl;
	else
	output << "precinctlistv01" <<endl;
	output << N << endl;
	output << "\tnb\tsp\tunshared\tarea\tpop\tvoteA\tvoteB\tcongD"<<endl;
	for (int i=0; i<N; i++){
		output << i;

		output << '\t';

		for (int j=0; j<pr[i].degree-1; j++){
			if (pr[i].neighbors[j]<N)
			output << pr[i].neighbors[j] << ',';
			else
			output << (-pr[i].self[j]-1) << ',';
		}
		if (pr[i].neighbors[pr[i].degree-1]<N)
		output << pr[i].neighbors[pr[i].degree-1];
		else
		output << (-pr[i].self[pr[i].degree-1]-1);

		output << '\t';

		for (int j=0; j<pr[i].degree-1; j++)
		output << pr[i].shared_perimeters[j] << ',';
		output << pr[i].shared_perimeters[pr[i].degree-1];
		output << '\t';
		output << 0;
		output << '\t';
		output << pr[i].area;
		output << '\t';
		output << pr[i].population;
		output << '\t';
		output << pr[i].voteA;
		output << '\t';
		output << pr[i].voteB;
		output << '\t';
		output << currentdistrict[i]+1;
		if (use_counties){
			output << '\t';
			output << pr[i].county;
		}
		output << endl;
	}
	output.close();
}





struct condbag{
	bool Perim;
	bool L1;
	bool L2;
	bool Popper;
	double perimthresh;
	double L1thresh;
	double L2thresh;
	double popperthresh;
	double popthresh;

	condbag(bool Pe,bool L,bool LL,bool Po, double pe, double l, double ll, double po, double pop){
		Perim=Pe;
		L1=L;
		L2=LL;
		Popper=Po;

		perimthresh=pe;
		L1thresh=l;
		L2thresh=ll;
		popperthresh=po;

		popthresh=pop; //population
	}
};



struct metric{
	double variance;
	double median_mean;
	double partisan_bias; // AMV EDIT #2: Added double partisan_bias;
	double BG_uniform; // AMV EDIT 10-26-19: Added uniform geometric bias
	double BG_modified; // AMV EDIT 11-1-19: Added modified geometric bias
	double seat_slide;
	double efficiency_gap;
	int seat_count;

	metric(double V, double M, double PB, double BGU, double BGM, double SS, double EG, int S){ // AMV EDIT #3: Added double PB, AMV EDIT 10-26-19 Added double BGU, AMV EDIT 11-1-19 Added double BGM
		variance=V;
		median_mean=M;
		partisan_bias=PB; // AMV EDIT #4: partisan_bias=PB;
		BG_uniform=BGU; // AMV EDIT 10-26-19: Added BG_uniform=BGU;
		BG_modified=BGM; // AMV EDIT 11-1-19: Added BG_modified=BGM;
		seat_slide=SS;
		efficiency_gap=EG;
		seat_count=S;
	}
	metric(){}
};

struct Ddata{
private:
	bool initialized;
public:
	int N;
	int * votesA;
	int * votesB;
	int * pops;
	double * Ashare;
	double * areas;
	double * perims;

	int * currentdistrict;
	rdpile<edge> edgeset;

	Ddata(){
		initialized=false;
	}
	

	Ddata(int NN){
		initialized=true;
		N=NN;
		votesA=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
		votesB=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
		Ashare= new double[g_NUMDISTRICTS];
		pops=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
		areas=new double[g_NUMDISTRICTS];   //district areas
		perims=new double[g_NUMDISTRICTS];  //district perimeters
		for (int k=0; k<g_NUMDISTRICTS; k++){
			votesA[k]=votesB[k]=pops[k]=0;
			areas[k]=perims[k]=0;
		}

		currentdistrict=new int[N+1];
		edgeset=rdpile<edge>(g_MAXEDGES,g_MAXEDGES*MAXDEGREE);        //set of edges on district boundaries
		
	}
	Ddata(const Ddata& source){ // This is a copy constructor. source is a Ddata object.
		if (source.initialized){
			initialized=true;
			N=source.N;
			votesA=new int[g_NUMDISTRICTS];      //running count of voteA in each district
			votesB=new int[g_NUMDISTRICTS];      //running count of voteB in each district
			Ashare= new double[g_NUMDISTRICTS];
			pops=new int[g_NUMDISTRICTS];      //running count of pop in each district
			areas=new double[g_NUMDISTRICTS];   //district areas
			perims=new double[g_NUMDISTRICTS];  //district perimeters
			
			for (int i=0; i<g_NUMDISTRICTS; i++){
				votesA[i]=source.votesA[i];
				votesB[i]=source.votesB[i];
				Ashare[i]=source.Ashare[i];
				pops[i]=source.pops[i];
				areas[i]=source.areas[i];
				perims[i]=source.perims[i];
			}
			currentdistrict=new int[N+1];
			for (int i=0; i<=N; i++)
			currentdistrict[i]=source.currentdistrict[i];

			edgeset=source.edgeset;
		}
		else
		initialized=false;
	}
	Ddata& operator= (const Ddata& source){
		if (&source!=this){
			if (source.initialized){
				initialized=true;
				N=source.N;
				votesA=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
				votesB=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
				Ashare= new double[g_NUMDISTRICTS];
				pops=new int[g_NUMDISTRICTS];      //running count of voteA-voteB in each district
				areas=new double[g_NUMDISTRICTS];   //district areas
				perims=new double[g_NUMDISTRICTS];  //district perimeters
				
				for (int i=0; i<g_NUMDISTRICTS; i++){
					votesA[i]=source.votesA[i];
					votesB[i]=source.votesB[i];
					Ashare[i]=source.Ashare[i];
					pops[i]=source.pops[i];
					areas[i]=source.areas[i];
					perims[i]=source.perims[i];
				}
				currentdistrict=new int[N+1];
				for (int i=0; i<=N; i++)
				currentdistrict[i]=source.currentdistrict[i];
				
				edgeset=source.edgeset;
			}
			else
			initialized=false;
		}
		return *this;
	}
	~Ddata(){
		if (initialized){
			delete[] votesA;
			delete[] votesB;
			delete[] Ashare;
			delete[] pops;
			delete[] areas;
			delete[] perims;
			delete[] currentdistrict;
		}
	}

};





struct datacount{
	int64_t variance_lessunusual;
	int64_t median_mean_lessunusual;
	int64_t partisan_bias_lessunusual; // AMV EDIT #5
	int64_t BG_uniform_lessunusual; // AMV EDIT 10-26-19
	int64_t BG_modified_lessunusual; // AMV EDIT 11-1-19
	int64_t seat_slide_lessunusual;
	int64_t efficiency_gap_lessunusual;
	int64_t seat_count_lessunusual;

	int64_t variance_moreunusual;
	int64_t median_mean_moreunusual;
	int64_t partisan_bias_moreunusual; // AMV EDIT #6
	int64_t BG_uniform_moreunusual; // AMV EDIT 10-26-19
	int64_t BG_modified_moreunusual; // AMV EDIT 11-1-19
	int64_t seat_slide_moreunusual;
	int64_t efficiency_gap_moreunusual;
	int64_t seat_count_moreunusual;

	int64_t totalsteps;

	datacount(int firstmore){
		variance_lessunusual=0;
		median_mean_lessunusual=0;
		partisan_bias_lessunusual=0; // AMV EDIT #7
		BG_uniform_lessunusual=0; // AMV EDIT 10-26-19
		BG_modified_lessunusual=0; // AMV EDIT 11-1-19
		seat_slide_lessunusual=0;
		efficiency_gap_lessunusual=0;
		seat_count_lessunusual=0;

		variance_moreunusual=firstmore;
		median_mean_moreunusual=firstmore;
		partisan_bias_moreunusual=firstmore; // AMV EDIT #8
		BG_uniform_moreunusual=firstmore; // AMV EDIT 10-26-19
		BG_modified_moreunusual=firstmore; // AMV EDIT 11-1-19
		seat_slide_moreunusual=firstmore;
		efficiency_gap_moreunusual=firstmore;
		seat_count_moreunusual=firstmore;

		totalsteps=firstmore;
	}

	datacount(){
		variance_lessunusual=0;
		median_mean_lessunusual=0;
		partisan_bias_lessunusual=0; // AMV EDIT #9
		BG_uniform_lessunusual=0; // AMV EDIT 10-26-19
		BG_modified_lessunusual=0; // AMV EDIT 11-1-19
		seat_slide_lessunusual=0;
		efficiency_gap_lessunusual=0;
		seat_count_lessunusual=0;

		variance_moreunusual=0;
		median_mean_moreunusual=0;
		partisan_bias_moreunusual=0; // AMV EDIT #10
		BG_uniform_moreunusual=0; // AMV EDIT 10-26-19
		BG_modified_moreunusual=0; // AMV EDIT 11-1-19
		seat_slide_moreunusual=0;
		efficiency_gap_moreunusual=0;
		seat_count_moreunusual=0;
		
		totalsteps=0;
	}
};



struct outputbag{
	bool doVariance;
	bool doMedianMean;
	bool doPartisanBias; // AMV EDIT #11
	bool doBGuniform; // AMV EDIT 10-26-19
	bool doBGmodified; // AMV EDIT 11-1-19
	bool doSeatSlide;
	bool doEfficiencyGap;
	bool doSeats;
	bool doHistogram;

	metric initial;

	int T;
	datacount * threadcounts;

	bool * thread_alive;

	int64_t ** reps_histogram;

public:


	outputbag(bool V, bool M, bool PB, bool BGU, bool BGM, bool SS, bool EG, bool S, bool H, double * Ashare, int * DvotesA, int* DvotesB, int initialvisits, int T){ // AMV EDIT #12: Added bool PB, AMV EDIT 10-26-19 Added bool BGU, AMV EDIT 11-1-19 Added bool BGM
		doVariance=V;
		doMedianMean=M;
		doPartisanBias=PB; // AMV EDIT #13
		doBGuniform=BGU; // AMV EDIT 10-26-19
		doBGmodified=BGM; // AMV EDIT 11-1-19
		doSeatSlide=SS;
		doEfficiencyGap=EG;
		doSeats=S;
		doHistogram=H;
		initial=metric(variance(Ashare), median_mean(Ashare), partisan_bias(Ashare), BG_uniform(Ashare), BG_modified(Ashare), seat_slide(Ashare), efficiency_gap(DvotesA,DvotesB), reps(Ashare)); // AMV EDIT #14: Added partisan_bias(Ashare), AMV EDIT 10-26-19 Added BG_uniform(Ashare), AMV EDIT 11-1-19 Added BG_modified(Ashare)

		this->T=T;
		threadcounts = new datacount[T];
		threadcounts[0]=datacount(initialvisits);
		for (int t=1; t<T; t++){
			threadcounts[t]=datacount();
		}

		thread_alive=new bool[T];
		for (int t=0; t<T; t++){
			thread_alive[t]=true;
		}
		
		reps_histogram = new int64_t*[T];
		for (int t=0; t<T; t++){
			reps_histogram[t]=new int64_t[g_NUMDISTRICTS+1];
			for (int i=0; i<=g_NUMDISTRICTS; i++)
			reps_histogram[t][i]=0;
		}

	}


	void initialoutput(){
		if (doVariance)
			cout << "initial variance is "<<initial.variance<<endl;
		if (doMedianMean)
			cout << "initial median_mean is "<<initial.median_mean<<endl;
		if (doPartisanBias)
			cout << "initial partisan_bias is " << initial.partisan_bias << endl;  // AMV EDIT #15
		if (doBGuniform)
			cout << "initial BG_uniform is " << initial.BG_uniform << endl; // AMV EDIT 10-26-19
		if (doBGmodified)
			cout << "initial BG_modified is " << initial.BG_modified << endl; // AMV EDIT 11-1-19
		if (doSeatSlide)
			cout << "initial seat_slide is "<<initial.seat_slide<<endl;
		if (doEfficiencyGap)
			cout << "initial efficiency_gap is "<<initial.efficiency_gap<<endl;
		if (doSeats)
			cout << "initial seat count is "<<initial.seat_count<<endl;
	}


	void update(int t, double* Ashare, int* DvotesA, int* DvotesB, int revisitations, int64_t STEPS){ //update counts for thread t
		assert(t<T);

		
		threadcounts[t].totalsteps+=revisitations;
		
		if (STEPS>0 && threadcounts[t].totalsteps>=STEPS && StatTest!=GCparallel){
			thread_alive[t]=false;               //thread will die
		}
		
		if (StatTest!=GCparallel){
			if (doVariance){
				if (variance(Ashare)>=initial.variance)
					threadcounts[t].variance_moreunusual+=revisitations;
				else
					threadcounts[t].variance_lessunusual+=revisitations;
			}
			if (doMedianMean){
				if (median_mean(Ashare)>=initial.median_mean) // This inequality direction is correct because the code calculates median_mean stupidly.
					threadcounts[t].median_mean_moreunusual+=revisitations;
				else
					threadcounts[t].median_mean_lessunusual+=revisitations;
			}
			 // AMV EDIT #16 ////
			if (doPartisanBias){
				if (partisan_bias(Ashare)<=initial.partisan_bias)
					threadcounts[t].partisan_bias_moreunusual+=revisitations;
				else
					threadcounts[t].partisan_bias_lessunusual+=revisitations;
			}
			/////////////////////
			// AMV EDIT 10-26-19 ////
			if (doBGuniform){
				if (BG_uniform(Ashare) >= initial.BG_uniform)
					threadcounts[t].BG_uniform_moreunusual+=revisitations;
				else
					threadcounts[t].BG_uniform_lessunusual+=revisitations;
			}
			/////////////////////////
			// AMV EDIT 11-1-19 ////
			if (doBGmodified){
				double BGMtest = BG_modified(Ashare);	// 11-23-21 edit: added flag for BG_modified - if an invalid seats-votes graph is produced (e.g. same Ashare in multiple districts), the measure won't get counted for this map
				if(BGMtest!=2){
					if (BG_modified(Ashare) >= initial.BG_modified)
						threadcounts[t].BG_modified_moreunusual+=revisitations;
					else
						threadcounts[t].BG_modified_lessunusual+=revisitations;
				}
			}
			/////////////////////////
			if (doSeatSlide){
				if (seat_slide(Ashare)<=initial.seat_slide)
				threadcounts[t].seat_slide_moreunusual+=revisitations;
				else
				threadcounts[t].seat_slide_lessunusual+=revisitations;
			}
			if (doEfficiencyGap){
				if (efficiency_gap(DvotesA,DvotesB)>=initial.efficiency_gap)
				threadcounts[t].efficiency_gap_moreunusual+=revisitations;
				else
				threadcounts[t].efficiency_gap_lessunusual+=revisitations;
			}
			if (doSeats){
				if (reps(Ashare)<=initial.seat_count)
				threadcounts[t].seat_count_moreunusual+=revisitations;
				else
				threadcounts[t].seat_count_lessunusual+=revisitations;
			}
			if (doHistogram){
				reps_histogram[t][reps(Ashare)]+=revisitations;
			}
		}
		else if  (StatTest==GCparallel && ((t>0 && threadcounts[t].totalsteps>=STEPS) || STEPS==0)){
			if (doVariance){
				if (variance(Ashare)>=initial.variance)
					threadcounts[t].variance_moreunusual+=1;
				else
					threadcounts[t].variance_lessunusual+=1;
			}
			if (doMedianMean){
				if (median_mean(Ashare)>=initial.median_mean)
					threadcounts[t].median_mean_moreunusual+=1;
				else
					threadcounts[t].median_mean_lessunusual+=1;
			}
			// AMV EDIT #17 //////////
			if (doPartisanBias){
				if (partisan_bias(Ashare)<=initial.partisan_bias)
					threadcounts[t].partisan_bias_moreunusual+=1;
				else
					threadcounts[t].partisan_bias_lessunusual+=1;
			}
			//////////////////////////
			// AMV EDIT 10-26-19 /////
			if (doBGuniform){
				if (BG_uniform(Ashare)>=initial.BG_uniform)
					threadcounts[t].BG_uniform_moreunusual+=1;
				else
					threadcounts[t].BG_uniform_lessunusual+=1;
			}
			//////////////////////////
			// AMV EDIT 11-1-19 ////
			if (doBGmodified){
				double BGMtest = BG_modified(Ashare);	// 11-23-21 edit: added flag for BG_modified - if an invalid seats-votes graph is produced (e.g. same Ashare in multiple districts), the measure won't get counted for this map
				if(BGMtest!=2){
					if (BG_modified(Ashare)>=initial.BG_modified)
						threadcounts[t].BG_modified_moreunusual+=1;
					else
						threadcounts[t].BG_modified_lessunusual+=1;
				}
			}
			//////////////////////////
			if (doSeatSlide){
				if (seat_slide(Ashare)<=initial.seat_slide)
					threadcounts[t].seat_slide_moreunusual+=1;
				else
					threadcounts[t].seat_slide_lessunusual+=1;
			}
			if (doEfficiencyGap){
				if (efficiency_gap(DvotesA,DvotesB)>=initial.efficiency_gap)
					threadcounts[t].efficiency_gap_moreunusual+=1;
				else
					threadcounts[t].efficiency_gap_lessunusual+=1;
			}
			if (doSeats){
				if (reps(Ashare)<=initial.seat_count){
					threadcounts[t].seat_count_moreunusual+=1;
				}
				else
					threadcounts[t].seat_count_lessunusual+=1;
			}
			if (doHistogram){
				reps_histogram[t][reps(Ashare)]+=1;
			}
		}
	}


	double pvalue(double epsilon){
		if (StatTest==rootep)
			return sqrt(2*epsilon);
			// AMV TO DO: Code a safeguard for when epsilon >= .5
		else if (StatTest==GCserial)
			return 2*epsilon;
		else if (StatTest==starsplit)
			return epsilon;
		else if (StatTest==GCparallel)
			return epsilon;
		else{
			cerr << "ERROR: unknown StatTest"<<endl;
			exit(-1);
		}
	}


	void  toscreen(int t, double * Ashare, int * DvotesA, int * DvotesB){
		assert(t<T);

		stringstream output;

		output << "Thread "<< t<< " reporting"<<endl;

		if (g_debuglevel>0){
			output << "Dshare is ";
			for (int j=0; j<g_NUMDISTRICTS; j++)
			output << setprecision(2) << 100*Ashare[j] <<"%, ";
			output << endl;
		}    


		int64_t Total_totalsteps=0;
		for (int t=0; t<T; t++)
		Total_totalsteps+=threadcounts[t].totalsteps;
		
		if (doHistogram){
			for (int j=0; j<=g_NUMDISTRICTS; j++){
				if (reps_histogram[t][j]>0){
					output << setw(2) << j << ": ";
					for (int k=0; k<(50*reps_histogram[t][j])/(threadcounts[t].totalsteps); k++)
					output << "*";
					for (int k=ceil(50*reps_histogram[t][j])/(threadcounts[t].totalsteps); k<50; k++)
					output << " ";
					output <<100*(double) reps_histogram[t][j]/(threadcounts[t].totalsteps)<<"%"<<endl;
				}
			}
			output << endl;
		}
		if (doVariance){
			output << "--FOR VARIANCE--"<<endl;
			output << "current variance is "<<setprecision(5)<< variance(Ashare)<<endl;
			output << "thread moreunusual is " << threadcounts[t].variance_moreunusual<<endl;
			output << "thread lessunusual is " << threadcounts[t].variance_lessunusual<<endl;

			int64_t Total_more=0;
			for (int t=0; t<T; t++)
			Total_more+=threadcounts[t].variance_moreunusual;

			double V_ep;
			if (StatTest==GCparallel){
				V_ep=((double) Total_more) / T;
				output << "total moreunusual/T is "<<Total_more<<"/"<<T<<endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				V_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is "<<Total_more<<"/"<<Total_totalsteps<<endl;
			}
			output << "ep="<<V_ep<<endl;
			output << "p="<<pvalue(V_ep)<<endl << endl;
		}
		if (doMedianMean){
			output << "--FOR MEDIAN/MEAN--"<<endl;
			output << "current median_mean is "<< setprecision(5)<< median_mean(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].median_mean_moreunusual<<endl;
			output << "thread lessunusual is " << threadcounts[t].median_mean_lessunusual<<endl;

			int64_t Total_more=0;
			for (int t=0; t<T; t++)
			Total_more+=threadcounts[t].median_mean_moreunusual;

			double M_ep;
			if (StatTest==GCparallel){
				M_ep=((double) Total_more) / T;
				output << "total moreunusual/T is "<<Total_more<<"/"<<T<<endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				M_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is "<<Total_more<<"/"<<Total_totalsteps<<endl;
			}

			output << "ep="<<M_ep<<endl;
			output << "p="<<pvalue(M_ep)<<endl << endl;
		}
		// AMV EDIT #18 ///////////////////////////////////////////////////////////////////////
		if (doPartisanBias){
			output << "--FOR PARTISAN_BIAS--" << endl;
			output << "current partisan_bias is " << setprecision(5) << partisan_bias(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].partisan_bias_moreunusual << endl;
			output << "thread lessunusual is " << threadcounts[t].partisan_bias_lessunusual << endl;
			
			int64_t Total_more=0;
			for (int t=0; t<T; t++)
				Total_more+=threadcounts[t].partisan_bias_moreunusual;
			
			double PB_ep;
			if (StatTest==GCparallel){
				PB_ep=((double) Total_more) / T;
				output << "total moreunusual/T is " << Total_more << "/" << T << endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				PB_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is " << Total_more << "/" << Total_totalsteps << endl;
			}
			
			output << "ep=" << PB_ep << endl;
			output << "p=" << pvalue(PB_ep) << endl << endl;
		}
		///////////////////////////////////////////////////////////////////////////////////////
		// AMV EDIT 10-26-19 ///////////////////////////////////////////////////////////////////////
		if (doBGuniform){
			output << "--FOR B_G UNIFORM--" << endl;
			output << "current BG_uniform is " << setprecision(5) << BG_uniform(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].BG_uniform_moreunusual << endl;
			output << "thread lessunusual is " << threadcounts[t].BG_uniform_lessunusual << endl;
			
			int64_t Total_more=0;
			for (int t=0; t<T; t++)
				Total_more+=threadcounts[t].BG_uniform_moreunusual;
			
			double BGU_ep;
			if (StatTest==GCparallel){
				BGU_ep=((double) Total_more) / T;
				output << "total moreunusual/T is " << Total_more << "/" << T << endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				BGU_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is " << Total_more << "/" << Total_totalsteps << endl;
			}
			
			output << "ep=" << BGU_ep << endl;
			output << "p=" << pvalue(BGU_ep) << endl << endl;
		}
		///////////////////////////////////////////////////////////////////////////////////////
		// AMV EDIT 11-1-19 ///////////////////////////////////////////////////////////////////////
		if (doBGmodified){
			output << "--FOR B_G MODIFIED--" << endl;
			output << "current BG_modified is " << setprecision(5) << BG_modified(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].BG_modified_moreunusual << endl;
			output << "thread lessunusual is " << threadcounts[t].BG_modified_lessunusual << endl;
			
			int64_t Total_more=0;
			for (int t=0; t<T; t++)
				Total_more+=threadcounts[t].BG_modified_moreunusual;
			
			double BGM_ep;
			if (StatTest==GCparallel){
				BGM_ep=((double) Total_more) / T;
				output << "total moreunusual/T is " << Total_more << "/" << T << endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				BGM_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is " << Total_more << "/" << Total_totalsteps << endl;
			}
			
			output << "ep=" << BGM_ep << endl;
			output << "p=" << pvalue(BGM_ep) << endl << endl;
		}
		///////////////////////////////////////////////////////////////////////////////////////
		
		if (doSeatSlide){
			output << "--FOR SEAT/SLIDE--"<<endl;
			output << "current seat_slide is "<< setprecision(5)<< seat_slide(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].seat_slide_moreunusual<<endl;
			output << "thread lessunusual is " << threadcounts[t].seat_slide_lessunusual<<endl;

			int64_t Total_more=0;
			for (int t=0; t<T; t++)
			Total_more+=threadcounts[t].seat_slide_moreunusual;

			double SS_ep;
			if (StatTest==GCparallel){
				SS_ep=((double) Total_more) / T;
				output << "total moreunusual/T is "<<Total_more<<"/"<<T<<endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				SS_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is "<<Total_more<<"/"<<Total_totalsteps<<endl;
			}
			output << "ep="<<SS_ep<<endl;
			output << "p="<<pvalue(SS_ep)<<endl << endl;
		}
		if (doEfficiencyGap){
			output << "--FOR Efficiency Gap--"<<endl;
			output << "current efficiency_gap is "<< setprecision(5)<< efficiency_gap(DvotesA,DvotesB) << endl;
			output << "thread moreunusual is " << threadcounts[t].efficiency_gap_moreunusual<<endl;
			output << "thread lessunusual is " << threadcounts[t].efficiency_gap_lessunusual<<endl;

			int64_t Total_more=0;
			for (int t=0; t<T; t++)
			Total_more+=threadcounts[t].efficiency_gap_moreunusual;

			double EG_ep;
			if (StatTest==GCparallel){
				EG_ep=((double) Total_more) / T;
				output << "total moreunusual/T is "<<Total_more<<"/"<<T<<endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				EG_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is "<<Total_more<<"/"<<Total_totalsteps<<endl;
			}      
			output << "ep="<<EG_ep<<endl;
			output << "p="<<pvalue(EG_ep)<<endl << endl;
		}
		if (doSeats){
			output << "--FOR Seat Count--"<<endl;
			output << "current seat_count is "<< setw(2)<< reps(Ashare) << endl;
			output << "thread moreunusual is " << threadcounts[t].seat_count_moreunusual<<endl;
			output << "thread lessunusual is " << threadcounts[t].seat_count_lessunusual<<endl;

			int64_t Total_more=0;
			for (int t=0; t<T; t++)
			Total_more+=threadcounts[t].seat_count_moreunusual;

			double S_ep;
			if (StatTest==GCparallel){
				S_ep=((double) Total_more) / T;
				output << "total moreunusual/T is "<<Total_more<<"/"<<T<<endl;
			}
			else{
				assert(StatTest==rootep || StatTest==GCserial || StatTest==starsplit);
				S_ep=((double) Total_more) / Total_totalsteps;
				output << "total moreunusual/steps is "<<Total_more<<"/"<<Total_totalsteps<<endl;
			}      

			output << "ep="<<S_ep<<endl;
			output << "p="<<pvalue(S_ep)<<endl << endl;
		}
		cout << output.str();
	}

	outputbag(const outputbag& source){
		doVariance=source.doVariance;
		doMedianMean=source.doMedianMean;
		doPartisanBias=source.doPartisanBias; // AMV EDIT #19
		doBGuniform=source.doBGuniform; // AMV EDIT 10-26-19
		doBGmodified=source.doBGmodified; // AMV EDIT 11-1-19
		doSeatSlide=source.doSeatSlide;
		doEfficiencyGap=source.doEfficiencyGap;
		doSeats=source.doSeats;
		doHistogram=source.doHistogram;
		initial=source.initial;
		T=source.T;
		
		threadcounts=new datacount[T];
		for (int t=0; t<T; t++){
			threadcounts[t]=source.threadcounts[t];
		}

		thread_alive=new bool[T];
		for (int t=0; t<T; t++){
			thread_alive[t]=source.thread_alive[t];
		}

		reps_histogram = new int64_t*[T];
		for (int t=0; t<T; t++){
			reps_histogram[t]=new int64_t[g_NUMDISTRICTS+1];
			for (int i=0; i<=g_NUMDISTRICTS; i++)
			reps_histogram[t][i]=source.reps_histogram[t][i];
		}
	}
	outputbag& operator= (const outputbag& source){
		if (&source!=this){
			doVariance=source.doVariance;
			doMedianMean=source.doMedianMean;
			doPartisanBias=source.doPartisanBias; // AMV EDIT #20
			doBGuniform=source.doBGuniform; // AMV EDIT 10-26-19
			doBGmodified=source.doBGmodified; // AMV EDIT 11-1-19
			doSeatSlide=source.doSeatSlide;
			doEfficiencyGap=source.doEfficiencyGap;
			doSeats=source.doSeats;
			doHistogram=source.doHistogram;
			initial=source.initial;
			T=source.T;
			
			threadcounts=new datacount[T];
			for (int t=0; t<T; t++){
				threadcounts[t]=source.threadcounts[t];
			}

			
			thread_alive=new bool[T];
			for (int t=0; t<T; t++){
				thread_alive[t]=source.thread_alive[t];
			}
			
			for (int t=0; t<T; t++){
				reps_histogram[t]=new int64_t[g_NUMDISTRICTS+1];
				for (int i=0; i<=g_NUMDISTRICTS; i++)
				reps_histogram[t][i]=source.reps_histogram[t][i];
			}
		}
		return *this;
	}
	~outputbag(){
		delete[] threadcounts;
		delete[] thread_alive;
		for (int t=0; t<T; t++)
		delete[] reps_histogram[t];
		delete[] reps_histogram;
	}

};




void chainstep(precinct * pr, edge e, Ddata &myD, double avgpop, condbag Test, bool Verify){  //add e.u(e.j) to district of e.u
	int Du=myD.currentdistrict[e.u];
	int v=pr[e.u].neighbors[e.j];
	int Dv=myD.currentdistrict[v];
	if (Du==Dv){
		cout <<"WTF! Du is "<<Du<<" and Dv is "<<Dv<<endl;
		exit(-1);
	}
	assert (Dv>=0); //not adding outside
	assert (Du!=Dv); 

	myD.votesA[Du]+= pr[v].voteA;	// Shift the precinct votes for Party A to new district
	myD.votesA[Dv]-= pr[v].voteA;
	myD.votesB[Du]+= pr[v].voteB;	// Shift the precinct votes for Party B to new district
	myD.votesB[Dv]-= pr[v].voteB;
	
	///////////////// 3-2-2020 EDITS /////////////////
	// Did we get a tie anywhere? If so, "toss a coin" to break the tie. 
	static mt19937_64 coin(93238462643383+g_seedoffset);
	static std::uniform_int_distribution<int> randflip(0,1);
	int tiebreaker_u=0;
	int voteshift_u=0;
	int tiebreaker_v=0;
	int voteshift_v=0;
	
	if (myD.votesA[Du]==myD.votesB[Du]){		// Check for tie in District u
		tiebreaker_u = randflip(coin);
		if (tiebreaker_u == 0){
			voteshift_u = -1;
		}
		else{	// tiebreaker == 1
			voteshift_u = 1;
		}
	}
	if (myD.votesA[Dv]==myD.votesB[Dv]){			// Check for tie in District v
		tiebreaker_v = randflip(coin);
		if (tiebreaker_v == 0){
			voteshift_v = -1;
		}
		else{	// tiebreaker == 1
			voteshift_v = 1;
		}
	}		// This ensures that no district will have a vote proportion of 0.5.
	
	myD.Ashare[Du]=((double) (myD.votesA[Du]+voteshift_u))/(myD.votesA[Du]+myD.votesB[Du]);	// Compute new Ashares
	myD.Ashare[Dv]=((double) (myD.votesA[Dv]+voteshift_v))/(myD.votesA[Dv]+myD.votesB[Dv]);

	/////////////// END 3-2-2020 EDITS ///////////////

	myD.pops[Du]+= pr[v].population;
	myD.pops[Dv]-= pr[v].population;

	myD.areas[Du]+= pr[v].area;
	myD.areas[Dv]-= pr[v].area;

	for (int l=0; l<pr[v].degree; l++){

		if (pr[v].neighbors[l]<0){  //outside
			myD.perims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
			myD.perims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv


		}
		else if (myD.currentdistrict[pr[v].neighbors[l]]==Du){
			if (!pr[pr[v].neighbors[l]].frozen){
				myD.edgeset.remove(edge(v,l));
				myD.edgeset.remove(edge(pr[v].neighbors[l],pr[v].self[l]));
			}
			myD.perims[Du]-=pr[v].shared_perimeters[l]; //because it intersects Du
			myD.perims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv
		}
		else if (myD.currentdistrict[pr[v].neighbors[l]]==Dv){
			if (!pr[pr[v].neighbors[l]].frozen){
				myD.edgeset.insert(edge(v,l));
				myD.edgeset.insert(edge(pr[v].neighbors[l],pr[v].self[l]));
			}
			myD.perims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
			myD.perims[Dv]+=pr[v].shared_perimeters[l]; //because it intersects Dv
		}
		else{ //avoids Du AND Dv
			myD.perims[Du]+=pr[v].shared_perimeters[l]; //because it avoids Du
			myD.perims[Dv]-=pr[v].shared_perimeters[l]; //because it avoids Dv
		}
	}

	myD.currentdistrict[v]=Du;


	if (Verify){     //CHECK whether we have to REVERSE the move
		bool LOOP=false;
		if (!validpop(myD.pops[Du],avgpop,Test.popthresh) || !validpop(myD.pops[Dv],avgpop,Test.popthresh)){
			LOOP=true;
			if (g_debuglevel>1)
			cout << "LOOP! for population violation"<<endl;
		}
		else if (Test.Perim && arraysum(myD.perims,g_NUMDISTRICTS)>Test.perimthresh){
			LOOP=true;
			if (g_debuglevel>1)
			cout << "LOOP! for Perimeter violation"<<endl;
		}
		else if (Test.Popper && (dcompact(myD.areas[Du],myD.perims[Du])>Test.popperthresh || dcompact(myD.areas[Dv],myD.perims[Dv])>Test.popperthresh)){
			LOOP=true;
			if (g_debuglevel>1)
			cout << "LOOP! for Polsby-Popper violation"<<endl;
		}
		else if (Test.L1 && compactsum(myD.areas,myD.perims,g_NUMDISTRICTS)>Test.L1thresh){
			LOOP=true;
			if (g_debuglevel>1)
			cout << "LOOP! for compactness L1 violation"<<endl;
		}
		else if (Test.L2 && compactL2(myD.areas,myD.perims,g_NUMDISTRICTS)>Test.L2thresh){
			LOOP=true;
			if (g_debuglevel>1)
			cout << "LOOP! for compactness L2 violation"<<endl;
		}

		//REVERSE the move if necessary
		if (LOOP){
			edge rev_e; //edge to reverse the move
			bool foundone=false;
			for (int l=0; l<pr[v].degree; l++){
				if (myD.currentdistrict[pr[v].neighbors[l]]==Dv){
					rev_e=edge(pr[v].neighbors[l],pr[v].self[l]);
					foundone=true;
					break;
				}
			}
			assert(foundone);
			chainstep(pr, rev_e, myD, avgpop, Test, false);
		}
	}
}


bool connectivity_check(precinct * pr, int * currentdistrict, rdpile<edge> & edgeset, edge e){

	int Du=currentdistrict[e.u];
	int v=pr[e.u].neighbors[e.j];
	int Dv=currentdistrict[v];
	
	int Dusegmentcount=0;
	int Dvsegmentcount=0;
	int precinct=pr[v].neighbors[0];

	int pointback=(pr[v].self[0]+1) % pr[precinct].degree;

	int lastincycle=pr[precinct].neighbors[pointback];
	int lastprecinct=lastincycle;
	int edgeindex=pr[precinct].self[pointback] ;
	
	
	int newprecinct=-1;
	int testcount=0;


	while ( ! ( newprecinct==pr[v].neighbors[0] && pr[lastprecinct].self[edgeindex]==pointback )  ){
		assert(precinct==pr[lastprecinct].neighbors[edgeindex]);
		testcount++;
		if (testcount>90)
			assert(testcount<100);
		
		int lastD=currentdistrict[lastprecinct];
		int newD=currentdistrict[precinct];

		if (lastD==Du && newD!=Du)
			Dusegmentcount++;
		else if (lastD==Dv && newD!=Dv)
			Dvsegmentcount++;

		edgeindex=(pr[lastprecinct].self[edgeindex]-1+pr[precinct].degree) % pr[precinct].degree;
		newprecinct=pr[precinct].neighbors[edgeindex];
		if (newprecinct==v){ //if I'm a neighbor of v, skip v for next step of cycle
			edgeindex=(edgeindex-1+pr[precinct].degree) % pr[precinct].degree;
			newprecinct=pr[precinct].neighbors[edgeindex];
		}
		lastprecinct=precinct;
		precinct=newprecinct;
	}

	
	assert(Dusegmentcount>0);
	if (Dvsegmentcount==0){
		cerr << "u is " <<e.u<<" and v is " <<v << endl;
	}
	assert(Dvsegmentcount>0);

	if (Dusegmentcount>1 || Dvsegmentcount>1)
	return false;
	else
	return true;
}



void MainLoop(Ddata &myD, outputbag & Out, condbag Tests, precinct * pr, int N, double avgpop,  gengetopt_args_info lineArgs, bool use_counties, int64_t STEPS, int t){
	mt19937_64 gen(31415926535897+t+g_seedoffset);
	int64_t period=pow(2,lineArgs.period_arg);
	int outputcount=0;
	int i=-1;
	while (Out.threadcounts[t].totalsteps<STEPS){
		i++;
		uniform_int_distribution<> intdist(0,myD.edgeset.count-1);
		int rindex=intdist(gen);
		edge e=myD.edgeset.access(rindex);    //we'll try adding vertex u(j) to u's district


		if (connectivity_check(pr, myD.currentdistrict, myD.edgeset, e)){  //we check connectivity before trying a move (other conds are CHECKed within chainstep after)
			chainstep(pr, e, myD, avgpop, Tests, true);
		}

		double p=((double) myD.edgeset.count)/g_MAXEDGES;
		std::geometric_distribution<int> revisits(p);
		int revisitations=1+revisits(gen);
		Out.update(t, myD.Ashare, myD.votesA, myD.votesB, revisitations, STEPS); // update counts for this step

		bool outputnow=false;
		if (! lineArgs.only_end_flag){
			//      if (i==0)
			//	outputnow=true;
			if (StatTest!=GCparallel &&  (Out.threadcounts[t].totalsteps/period > (Out.threadcounts[t].totalsteps-revisitations)/period)){
				outputnow=true;
				for (int s=0; s<t; s++){
					if (Out.thread_alive[s]==true)
					outputnow=false;  //output only if we are smallest alive thread 
				}
			}
			if (Out.threadcounts[t].totalsteps>=STEPS && StatTest!=GCparallel)
			outputnow=true;
		}
		if (outputnow){
			outputcount++;
			Out.toscreen(t, myD.Ashare, myD.votesA, myD.votesB);
			
			if (lineArgs.stages_flag){
				if (lineArgs.svg_filename_given){
					char svgfilename[100];
					sprintf(svgfilename,"%s_%d_%d.svg",lineArgs.svg_filename_arg,t,outputcount);
					char inputsvgfilename[100];
					sprintf(inputsvgfilename,"%s",lineArgs.inputsvg_filename_arg);
					tosvg(svgfilename,inputsvgfilename,pr,myD.currentdistrict,lineArgs.svg_firstline_arg,N);
				}
				
				if (lineArgs.precinct_filename_given){
					char precinctfilename[100];
					sprintf(precinctfilename,"%s_%d_%d",lineArgs.precinct_filename_arg,t,outputcount);
					tofile(precinctfilename,pr,myD.currentdistrict,N,use_counties);
				}
			}
		}
	}
	if (lineArgs.svg_filename_given){
		char svgfilename[100];
		sprintf(svgfilename,"%s_%d_%d.svg",lineArgs.svg_filename_arg,t,outputcount);
		char inputsvgfilename[100];
		sprintf(inputsvgfilename,"%s",lineArgs.inputsvg_filename_arg);
		tosvg(svgfilename,inputsvgfilename,pr,myD.currentdistrict,lineArgs.svg_firstline_arg,N);
	}
	
	if (lineArgs.precinct_filename_given){
		char precinctfilename[100];
		sprintf(precinctfilename,"%s_%d_%d",lineArgs.precinct_filename_arg,t,outputcount);
		tofile(precinctfilename,pr,myD.currentdistrict,N,use_counties);
	}
}




int main(int argc, char* argv[])
{
	gengetopt_args_info lineArgs;
	if (cmdline_parser(argc, argv, &lineArgs)) {
		cmdline_parser_print_help();
		return 1;
	}

	g_seedoffset=lineArgs.seed_arg;	// default zero

	int T=1; //threads

	if (!strcmp(lineArgs.statistic_arg,"rootep"))   {
		StatTest=rootep;
	}
	else if (!strcmp(lineArgs.statistic_arg,"GCserial"))   {
		StatTest=GCserial;
		T=2;
	}
	else if (!strcmp(lineArgs.statistic_arg,"starsplit"))   {
		StatTest=starsplit;
		T=lineArgs.threads_arg;
	}
	else if (!strcmp(lineArgs.statistic_arg,"GCparallel"))   {
		StatTest=GCparallel;
		if  (!lineArgs.branches_given){
			cerr << "ERROR: --branches required for GCparallel"<<endl;
			exit(-1);
		}
		T=pow(2,lineArgs.branches_arg);
	}
	else{
		cerr << "Unknown Statistical Test requested!!\n";
		exit(-1);
	}

	g_NUMDISTRICTS=lineArgs.numdists_arg;

	int numberfrozen=lineArgs.freeze_given;
	bool * frozen_districts;  //which districts are frozen
	frozen_districts=new bool[g_NUMDISTRICTS];

	for (int k=0; k<g_NUMDISTRICTS; k++)
	frozen_districts[k]=0;                       //zero them

	for (int i=0; i<numberfrozen; i++)
	frozen_districts[lineArgs.freeze_arg[i]-1]=1;//set to one the given frozen 
	//(indexing from 0)

	condbag Tests(lineArgs.perimeter_given, lineArgs.L1_compactness_given, lineArgs.L2_compactness_given, lineArgs.polsby_popper_given, lineArgs.perimeter_arg, lineArgs.L1_compactness_arg, lineArgs.L2_compactness_arg, lineArgs.polsby_popper_arg, lineArgs.poperror_arg);


	g_debuglevel=0;

	int N; //num precincts
	precinct * pr; //array of precincts


	ifstream myfile (lineArgs.filename_arg);
	string line;
	if (!myfile.good()){
		cerr << "ERROR with file"<<endl;
		exit(-1);
	}
	bool use_counties=false;
	string old_header ("precinctlistv01");
	string county_header ("precinctlistv02");
	getline (myfile,line); //file header
	if (old_header.compare(line)){
		if (county_header.compare(line)){
			cerr << "ERROR: incorrect file header at pos " << county_header.compare(line) << endl;
			exit(-1);
		}
		use_counties=true;
	}
	if (lineArgs.counties_flag && ! use_counties){
		cerr << "ERROR: preserving preserved counties required a file v02 file with county data!"<<endl;
		exit(-1);
	}
	getline(myfile,line); //#precints
	N=atoi(line.c_str());
	pr=new precinct[N+1];  


	if (! lineArgs.only_end_flag)
	cout << "We have "<<N<<" precincts."<<endl;			// N = Number of Precincts
	//upper bound on crossing edges.					// What are "crossing edges"?
	g_MAXEDGES=2*(2*N+(g_NUMDISTRICTS-numberfrozen)-6);	// Why subtract 6?!

	getline(myfile,line);  //skip data header

	int I=0;
	int sum=0;
	while (getline(myfile,line) && I<N){
		vector<string> dataline;
		Tokenize(line.c_str(), dataline, "\t", true);
		if (dataline.size()!=9+use_counties){
			cerr << "ERROR: shouldn't there be "<<9+use_counties<<" chunks per line?"<<endl;
			cerr << "there are "<<dataline.size()<<" at I="<<I<<endl;
			exit(-1);
		}
		if (atoi(dataline[0].c_str())!=I){
			cerr << atoi(dataline[0].c_str()) << endl;
			cerr << I << endl;
			cerr << "ERROR: line numbers don't match"<<endl;
			exit(-1);
		}
		pr[I]=precinct(dataline,lineArgs.flip_flag,use_counties);

		sum+=pr[I].degree;
		I++;
	}

	if (I<N || !myfile.eof() ){
		cerr << "ERROR: Line numbers don't match precinct count"<<endl;
		exit(-1);
	}

	if (lineArgs.counties_flag && use_counties)
	freezedistrictsbycounty(pr,N);

	//create backneighbor lists:
	computeselves(pr,N);



	int ** adjM;                                       //precinct adjacency matrix

	Ddata firstD(N);
	mt19937_64 gen(314159265358979+g_seedoffset);

	for (int i=0; i<=N; i++)
	firstD.currentdistrict[i]=pr[i].original_district;


	adjM=new int*[N];
	for (int i=0; i<N; i++){
		adjM[i]=new int[N];
		for (int j=0; j<N; j++)
		adjM[i][j]=0;
	}

	for (int i=0; i<N; i++){
		for (int j=0; j<pr[i].degree; j++){
			int idistrict=pr[i].original_district;
			int jdistrict=pr[pr[i].neighbors[j]].original_district;
			if (pr[i].neighbors[j]<N){  //not state boundary connection
				adjM[i][pr[i].neighbors[j]]++;
				if (idistrict!=jdistrict &&
						frozen_districts[idistrict]==0 &&
						frozen_districts[jdistrict]==0 &&
						!pr[i].frozen &&
						!pr[pr[i].neighbors[j]].frozen
						){
					firstD.edgeset.insert(edge(i,j)); 
				}
			}
		}
	}
	if (! lineArgs.only_end_flag)
	cout << "There are "<<firstD.edgeset.count<<"/2 boundary edges"<<endl;
	assert(firstD.edgeset.count%2==0);



	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			if (adjM[i][j]!=adjM[j][i]){
				cerr << "ERROR: graph is not undirected!"<< endl;
				cerr << "bad pair is "<<i<<","<<j<<"."<<endl;
				exit(-1);
			}
		}
	}


	/////////////INITIAL CONNECTIVITY CHECK///////
	int count=0;
	if (g_debuglevel>1)
	cout << "connectivity check..."<<endl;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		int localcount=0;
		queue <int> BFSqueue;
		for (int i=0; i<N; i++){
			if (pr[i].original_district==k){
				BFSqueue.push(i);
				pr[i].giant=true;
				count++; localcount++;
				break;
			}
		}
		while (!BFSqueue.empty()){
			int myindex=BFSqueue.front();
			BFSqueue.pop();
			assert(pr[myindex].giant==true);
			for (int j=0; j<pr[myindex].degree; j++){
				if (pr[pr[myindex].neighbors[j]].original_district==pr[myindex].original_district && pr[pr[myindex].neighbors[j]].giant==false){
					pr[pr[myindex].neighbors[j]].giant=true; count++; localcount++;
					BFSqueue.push(pr[myindex].neighbors[j]);
				}
			}
		}
		if (g_debuglevel>1)
		cerr << "district "<<k<< " giant has "<<localcount<<" precincts."<<endl;
	}

	bool passed=true;
	for (int i=0; i<N; i++){
		if (pr[i].giant==false){
			passed=false;
			cerr << "ERROR: Precinct "<<i<< " is not in its giant!"<<endl;
		}
	}

	if (passed){
		assert(count==N);
		if (g_debuglevel>1)
		cout << "passed connectivity check" <<endl;
	}
	else
	exit(-1);
	///////////END INITIAL CONNECTIVITY CHECK//////////////

	int64_t steps=pow(2,lineArgs.steps_arg);




	//filling in first step...

	int revisitations;
	{
		double p=((double) firstD.edgeset.count)/g_MAXEDGES;
		std::geometric_distribution<int> revisits(p);
		revisitations=1+revisits(gen);
		if (g_debuglevel>0){
			cout << "first edgeset.count is "<<firstD.edgeset.count<<endl;
		}
		if (! lineArgs.only_end_flag){
			cout << "first p is "<<p<<endl;
			cout << "first revisitations is "<<revisitations<<endl;
		}
	}
	int totalpop=0;
	for (int i=0; i<N; i++){
		firstD.pops[firstD.currentdistrict[i]]+=pr[i].population;
		totalpop+=pr[i].population;
		firstD.votesA[firstD.currentdistrict[i]]+= pr[i].voteA;
		firstD.votesB[firstD.currentdistrict[i]]+= pr[i].voteB;
		firstD.areas[firstD.currentdistrict[i]]+=  pr[i].area;
		for (int j=0; j<pr[i].degree; j++){
			if (pr[i].neighbors[j]>=0 && firstD.currentdistrict[i]!=firstD.currentdistrict[pr[i].neighbors[j]]){
				firstD.perims[firstD.currentdistrict[i]]+=pr[i].shared_perimeters[j];
			}
			else if (pr[i].neighbors[j]<0){
				firstD.perims[firstD.currentdistrict[i]]+=pr[i].shared_perimeters[j];
			}
		}
	}
	double avgpop = (double) totalpop/g_NUMDISTRICTS;
	int Avotes=0;
	int Bvotes=0;
	for (int k=0; k<g_NUMDISTRICTS; k++){
		assert(firstD.votesA[k]!=firstD.votesB[k]);
		firstD.Ashare[k]=((double) firstD.votesA[k])/(firstD.votesA[k]+firstD.votesB[k]);
		Avotes+=firstD.votesA[k];
		Bvotes+=firstD.votesB[k];
		if (g_debuglevel>0)
		cout << "district "<<k<<" has " <<firstD.votesA[k]<<" voteA and "<<firstD.votesB[k]<<" voteB.  share is "<<firstD.Ashare[k]<<endl;
	}
	if (! lineArgs.only_end_flag){
		cout << "A has "<<Avotes <<" votes"<<endl;
		cout << "B has "<<Bvotes <<" votes"<<endl;
		cout << "total population is "<<totalpop<<endl;
	}
	 // AMV EDIT #21: Added lineArgs.partisan_bias_flag,
	 // AMV EDIT 10-26-19: Added lineArgs.BG_uniform_flag,
	 // AMV EDIT 11-1-19: Added lineArgs.BG_modified_flag,
	outputbag Out(lineArgs.variance_flag, lineArgs.median_mean_flag, lineArgs.partisan_bias_flag, lineArgs.BG_uniform_flag, lineArgs.BG_modified_flag, lineArgs.seat_slide_flag, lineArgs.efficiency_gap_flag, lineArgs.seats_flag,lineArgs.histogram_flag, firstD.Ashare, firstD.votesA, firstD.votesB, revisitations, T);

	if (! lineArgs.only_end_flag)
	Out.initialoutput();


	if (Tests.Perim){
		for (int i=0; i<N; i++)
		if (arraysum(firstD.perims,g_NUMDISTRICTS)>Tests.perimthresh){
			cerr << "ERROR: initial districting violates perimeter criterion.  Current perimeter is "<< arraysum(firstD.perims,g_NUMDISTRICTS)  <<endl;
			exit(-1);
		}
	}
	if (Tests.Popper){
		for (int k=0; k<g_NUMDISTRICTS; k++){
			if (dcompact(firstD.areas[k],firstD.perims[k])>Tests.popperthresh){
				cerr << "ERROR: given district "<<k<<" does not satisfy compactness requirement!\n It has compactness value "<< dcompact(firstD.areas[k],firstD.perims[k]) <<"."<<endl;
				exit(-1);
			}
		}
	}
	if (Tests.L1){
		for (int i=0; i<N; i++)
		if (compactsum(firstD.areas,firstD.perims,g_NUMDISTRICTS)>Tests.L1thresh){
			cerr << "ERROR: initial districting violates L1 compactness criterion.  Current L1 norm is "<< compactsum(firstD.areas,firstD.perims,g_NUMDISTRICTS)  <<endl;
			exit(-1);
		}
	}
	if (Tests.L2){
		for (int i=0; i<N; i++)
		if (compactL2(firstD.areas,firstD.perims,g_NUMDISTRICTS)>Tests.L2thresh){
			cerr << "ERROR: initial districting violates L2 criterion.  Current L2 norm is "<< compactL2(firstD.areas,firstD.perims,g_NUMDISTRICTS)  <<endl;
			exit(-1);
		}
	}



	//MAIN LOOP
	//  double p=((double) firstD.edgeset.count)/g_MAXEDGES;
	if (StatTest==rootep){
		int t=0;
		Ddata myD=firstD;
		MainLoop(myD, Out, Tests, pr, N, avgpop, lineArgs, use_counties, steps, t);
		Out.toscreen(0, myD.Ashare, myD.votesA, myD.votesB);
	}
	else if (StatTest==GCserial){
#pragma omp parallel for
		for (int t=0; t<2; t++){
			Ddata myD=firstD;
			MainLoop(myD, Out, Tests, pr, N, avgpop, lineArgs, use_counties, steps, t);
			Out.toscreen(0, myD.Ashare, myD.votesA, myD.votesB);
		}
	}
	else if (StatTest==starsplit){
		uniform_int_distribution<int64_t> Lintdist(0,steps-1);
		int64_t * threadsteps;
		threadsteps=new int64_t[T];
		for (int t=0; t<T; t++)
		threadsteps[t]=steps;
		threadsteps[0]=Lintdist(gen);
		threadsteps[1]=steps-threadsteps[0];
		if (! lineArgs.only_end_flag){
			cout << "We'll branch after "<<threadsteps[0]<<" steps"<<endl;
		}
		Ddata * threadD;
		threadD=new Ddata[T];
		threadD[0]=firstD;    //First threads begin from initial state;
		threadD[1]=firstD;
#pragma omp parallel for
		for (int t=0; t<2; t++){
			MainLoop(threadD[t], Out, Tests, pr, N, avgpop, lineArgs, use_counties, threadsteps[t], t);
		}
		for (int t=2; t<T; t++){
			threadD[t]=threadD[0]; //remaining threads branch from end of one of the previous threads
		}
#pragma omp parallel for
		for (int t=2; t<T; t++){
			MainLoop(threadD[t], Out, Tests, pr, N, avgpop, lineArgs, use_counties, threadsteps[t], t);
		}
		Out.toscreen(T-1, threadD[T-1].Ashare, threadD[T-1].votesA, threadD[T-1].votesB);
	}
	else if (StatTest==GCparallel){
		int t=0;
		Ddata * threadD=new Ddata[T];
		threadD[0]=firstD;    //First threads begin from initial state;
		MainLoop(threadD[0], Out, Tests, pr, N, avgpop, lineArgs, use_counties, steps, t);
		Ddata outD;
#pragma omp parallel for
		for (int t=1; t<T; t++){
			threadD[t]=threadD[0];
			MainLoop(threadD[t], Out, Tests, pr, N, avgpop, lineArgs, use_counties, steps, t);
			if (!lineArgs.only_end_flag)
				Out.toscreen(t, threadD[t].Ashare, threadD[t].votesA, threadD[t].votesB);
		}
		Out.toscreen(1, threadD[1].Ashare, threadD[1].votesA, threadD[1].votesB);
	}

}


/*  if (lineArgs.svg_filename_given){
	char svgfilename[100];
	sprintf(svgfilename,"%s.svg",lineArgs.svg_filename_arg);
	char inputsvgfilename[100];
	sprintf(inputsvgfilename,"%s",lineArgs.inputsvg_filename_arg);
	tosvg(svgfilename,inputsvgfilename,pr,currentdistrict,lineArgs.svg_firstline_arg,N);
}

if (lineArgs.precinct_filename_given){
	char precinctfilename[100];
	sprintf(precinctfilename,"%s",lineArgs.precinct_filename_arg);
	tofile(precinctfilename,pr,currentdistrict,N,use_counties);
	}*/
