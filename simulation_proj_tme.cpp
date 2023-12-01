#include <iostream>
#include <random>
// #include <opencv2/opencv.hpp>
#include <time.h>

using namespace std;

const int NULL_CELL_TYPE = 0; //signifies an absence of cell
const int BLOOD_VESSEL_TYPE = 1;
const int CANCER_STEM_CELL_TYPE = 2;
const int NORMAL_CELL_TYPE = 3;
const int CANCER_CELL_TYPE = 4;

//probabilties of celltypes for initialization
const double prob_cancer_cell = 0.15;
const double prob_cancer_stem_cell=0.05;
const double prob_normal_cell =0.80;

const double prob_stem_cell_divides_to_cancer_cell = 0.3;

//indices for chemical entity array
const int VEGF_INDEX = 0;
const int O2_INDEX = 1;
const int GLUCOSE_INDEX =2;
const int LACTIC_INDEX = 3;


const double TUMOR_VEGF_CONC = 25.0; //measured in ppm ref:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3672077/
const double NORMAL_VEGF_CONC = 3.0; //measured in ppm
const double BLOOD_VESSEL_O2_CONC = 30.0; //measured in mmhg partial pressure
const double TUMOUR_PH_CONC = 6.6; //Mesured in PH
const double NORMAL_PH_CONC = 7.0;

const double NORMAL_LACTIC_ACID_CONC = 0.0;// measured in ppm
const double TUMOR_LACTIC_ACID_CONC = 7.0;

const double BLOOD_VESSEL_GLUCOSE_CONC = 100; //TODO - lets assume it is measured in PPM
const int BLOOD_VESSEL_BASE = 2; // number of rows in cell type matrix which pose as the major artery (starting from 0th row.)


//Cell Sizes 
//reference: https://www.pnas.org/doi/10.1073/pnas.2303077120
//Mass estimate in grams
const double SIZE_MACROPHAGE = 5.5e-8;
const double SIZE_EPITHELIAL_CELL = 8.0e-9;
const double SIZE_T_LYMPHOCYTE = 1.25e-10;
const double SIZE_DENDRITIC_CELL = 5.5e-9;
const double SIZE_RBC_CELL = 9.79e-11; //RBCs are 7.5e-6 in dia as per https://publiclab.org/notes/MaggPi/08-19-2019/microscope-calibration-with-image-sequencer


//// Global Random Number Generators
default_random_engine generator;


//// When a cell takes an action the follwing CellAction object is returned.
//// The calling fucntion can update the global state and environment and
/// list of cells accordingly.

const int CELL_ACTION_NULL =0;
const int CELL_ACTION_DIE = 1;
const int CELL_ACTION_DIVIDE = 2;
const int CELL_ACTION_MOVE = 3;
struct CellAction
{
	int action_type; 
	int xpos;
	int ypos;
	int cell_type;
	CellAction(int atyp=CELL_ACTION_NULL, int x=0, int y=0, int ctyp=NULL_CELL_TYPE)
	{
		action_type = atyp;
		xpos = x;
		ypos = y;
		cell_type = ctyp;

	}
	CellAction& operator = (const CellAction& rhs)
	{
		action_type = rhs.action_type;
		xpos = rhs.xpos;
		ypos = rhs.ypos;
		cell_type = rhs.cell_type;

		return *this;
	}
};

class Cell
{
protected:
	int ndivisions;

	// the weights of action is updated based on sesing the environment
	// then with epsilon probability one of actions is chosen and (1-epsilon prob) no action is done
	double prob_wt_div;
	double prob_wt_move;
	double prob_wt_die;
	double epsilon;

	int sensing_window;

	//values sensed over time
	//Use same chem entity index for the array 
	vector<vector<double> > sensed_chem_vals;

	double avg_val(double* conc_array, const int w, const int h)
	{
		double avg = 0.0;
		int n =1;
		avg = conc_array[y*w +x];
		if (x > 1)
		{
			avg += conc_array[y*w +x-1];
			n +=1;
		}
		if ( x< w-1)
		{
			avg += conc_array[y*w +x+1];
			n +=1;
		}
		if (y>1)
		{
			avg += conc_array[(y-1)*w +x];
			n +=1;
		}
		if (y < h-1)
		{
			avg += conc_array[(y+1)*w +x];
			n +=1;
		}
		avg /= n;
		return avg;
	}

	// double prob_die(const int timestep)
	// {
	// 	const double a1 = 1.0;
	// 	const double a2 = 0.5;
	// 	const double a3 = 1.0;
	// 	const double b = 1.0;

	// 	double avgo2 = 0.0;	
	// 	for (int i=0; i< sensed_chem_vals[O2_INDEX].size(); i++)
	// 	{
	// 		avgo2 += sensed_chem_vals[O2_INDEX][i];
	// 	}
	// 	avgo2 /= sensed_chem_vals[O2_INDEX].size();

	// 	double avg_glucose = 0.0;
	// 	for (int i=0; i< sensed_chem_vals[GLUCOSE_INDEX].size(); i++)
	// 	{
	// 		avg_glucose += sensed_chem_vals[GLUCOSE_INDEX][i];
	// 	}
	// 	avg_glucose /= sensed_chem_vals[GLUCOSE_INDEX].size();


	// 	const double pd = 1.0/(1.0+b*exp(-a1*(avgo2) -a2*(avg_glucose) +a3*(timestep)));
	// 	return pd;
	// }

public:

	int x;
	int y;
	int cell_type;

	bool just_born;
	bool dead_cell;

	Cell(const int xloc, const int yloc, const int typ,
		unsigned char* cell_type_map, const int w, const int h, const int n_chem_entities)
	{
		sensing_window = 5;
		x = xloc;
		y = yloc;
		cell_type = typ;
		ndivisions =0;
		cell_type_map[yloc*w +xloc] =cell_type;
		just_born = true;
		dead_cell = false;

		epsilon = 1.0; //No action by default

		for (int c =0; c <n_chem_entities; c++)
		{
			sensed_chem_vals.push_back(vector<double>());
		}
		
	}

	void Sense(const vector<double*>& conc_array,
				const int w, const int h,
				const int timestep)
	{
		// Sense the environment
		double avg_o2 =0.0;
		double avg_lactate = 0.0;
		double avg_glucose = 0.0;
		double avg_vegf = 0.0;

		avg_o2 = avg_val(conc_array[O2_INDEX],w,h);
		sensed_chem_vals[O2_INDEX].push_back(avg_o2);
		if (sensed_chem_vals[O2_INDEX].size() > sensing_window)
		{
			sensed_chem_vals[O2_INDEX].erase(sensed_chem_vals[O2_INDEX].begin());
		}

		avg_lactate = avg_val(conc_array[LACTIC_INDEX],w,h);
		sensed_chem_vals[LACTIC_INDEX].push_back(avg_lactate);
		if (sensed_chem_vals[LACTIC_INDEX].size() > sensing_window)
		{
			sensed_chem_vals[LACTIC_INDEX].erase(sensed_chem_vals[LACTIC_INDEX].begin());
		}

		avg_glucose = avg_val(conc_array[GLUCOSE_INDEX],w,h);
		sensed_chem_vals[GLUCOSE_INDEX].push_back(avg_glucose);
		if (sensed_chem_vals[GLUCOSE_INDEX].size() > sensing_window)
		{
			sensed_chem_vals[GLUCOSE_INDEX].erase(sensed_chem_vals[GLUCOSE_INDEX].begin());
		}

		avg_vegf = avg_val(conc_array[VEGF_INDEX],w,h);
		sensed_chem_vals[VEGF_INDEX].push_back(avg_vegf);
		if (sensed_chem_vals[VEGF_INDEX].size() > sensing_window)
		{
			sensed_chem_vals[VEGF_INDEX].erase(sensed_chem_vals[VEGF_INDEX].begin());
		}

	}

	//returns the action done - so that if it is die, then that cell can be erased
	CellAction ChooseActionAndAct(vector<Cell>& cell_vec, 
				unsigned char* cell_type_map, 
				const int num_chem_entities,
				const vector<double*>& conc_array,
				const int w, const int h,
				const int timestep)
	{
	
	// 	// the weights of action is updated based on sesing the environment
	// // then with (1-epsilon) probability one of actions is chosen and epsilon prob no action is done
	// 	double prob_wt_div;
	// 	double prob_wt_move;
	// 	double prob_wt_die;
	// 	double epsilon;

		CellAction d;
		d.action_type = CELL_ACTION_NULL;

		uniform_real_distribution<double> dist(0.0,1.0);
		double p = dist(generator);

		//cout <<"p =" <<p <<" epsilon= "<<epsilon <<"\n";

		if (p <epsilon)
		{
			//No action
			//cout <<"NO action: \n";
			d.action_type = CELL_ACTION_NULL;
			return d;
		}
		vector<double> act_prob_wt_vec;
		act_prob_wt_vec.push_back(prob_wt_div);
		act_prob_wt_vec.push_back(prob_wt_move);
		act_prob_wt_vec.push_back(prob_wt_die);
		discrete_distribution<size_t> action_type_distr(act_prob_wt_vec.begin(), act_prob_wt_vec.end());
		const int a = action_type_distr(generator);

		//cout <<"action: " <<a <<"\n";
		
		if (a ==0)
		{
			d = Divide(cell_vec,cell_type_map, w, h,num_chem_entities);
		}
		else if (a == 1)
		{
			d = Move(cell_vec,cell_type_map, w, h,num_chem_entities);
		}
		else if ( a == 2)
		{
			d = Die(cell_vec,cell_type_map, w, h,num_chem_entities);
		}
		else
		{
			cout <<" Bad value for cell action \n";
		}

		return d;

	}

	CellAction Die(vector<Cell>& cell_vec,unsigned char* cell_type_map, const int w, const int h,
		const int num_chem_entities)
	{
		CellAction d(CELL_ACTION_DIE,0,0,cell_type);
		cell_type_map[y*w+x] =NULL_CELL_TYPE;

		dead_cell = true;
		return d;
	}

	pair<int,int> FindPosToMove(unsigned char* cell_type_map, const int w, const int h)
	{
		//randomly choose an empty spot in the neighborhood. Neighborhood size dependent on 
		// lactic acid levels.

		int ynew = 0;
		int xnew =0; 

		double avg_lactic = 0.0;
		if (sensed_chem_vals[LACTIC_INDEX].size() >0)
		{
			for (int i=0; i< sensed_chem_vals[LACTIC_INDEX].size(); i++)
			{
				avg_lactic += sensed_chem_vals[LACTIC_INDEX][i];
			}
			avg_lactic /= sensed_chem_vals[LACTIC_INDEX].size();
		}

		//cout << "avg_lactic = " << avg_lactic <<" NORMAL_LACTIC_ACID_CONC = " <<NORMAL_LACTIC_ACID_CONC <<"\n";
		
		int nd = 2 + (avg_lactic - NORMAL_LACTIC_ACID_CONC)/3;

		int xs = x-nd;
		if (xs <1)
			xs =1;
		int xe = x +nd;
		if (xe >= w-1)
			xe = w-1;

		int ys = y-nd;
		if (ys <1)
			ys =1;
		int ye = y +nd;
		if (ye >= h-1)
			ye = h-1;

		uniform_int_distribution<> x_pos_dist(x-nd,x+nd);
		uniform_int_distribution<> y_pos_dist(y-nd,y+nd);

		bool found_pos = false;
		int niter =0;
		int max_iters = 25+ (2*nd+1)*(2*nd +1);
		//cout <<"nd = " <<nd <<"\n";
		do {
			xnew = x_pos_dist(generator);
			ynew = y_pos_dist(generator);

			if (cell_type_map[ynew*w +xnew] == NULL_CELL_TYPE)
			{
				found_pos =true;
			}
			niter++;

			//cout <<"niter = " <<niter <<"\n";

		}while(!found_pos && (niter <max_iters));

		if (cell_type_map[ynew*w +xnew] != NULL_CELL_TYPE)
		{
			ynew = -1;
			xnew = -1;
		}
		return make_pair(xnew,ynew);
		

	}

	CellAction  Divide(vector<Cell>& cell_vec,unsigned char* cell_type_map, const int w, const int h, 
				const int num_chem_entities)
	{
		CellAction d(CELL_ACTION_DIVIDE,0,0,0);

		/// TODO - Logic on where new cell  will move still not written

		pair<int,int> newpos = FindPosToMove(cell_type_map,w, h);
		int xnew = newpos.first;
		int ynew = newpos.second;

		if ( (xnew < 1) || (ynew <1) ||
			 (xnew >= w-1) ||  (ynew >= h-1) )
		{
			d.action_type=CELL_ACTION_NULL;
			d.xpos = -1;
			d.ypos = -1;	
			return d;
		}

		d.xpos = xnew;
		d.ypos = ynew;

		d.cell_type = cell_type;
		//new cell type can vary if the type is a cancer stem cell
		if (cell_type == CANCER_STEM_CELL_TYPE)
		{
			//with some probailty we get a new cancer stem cell
			uniform_real_distribution<double> dist(0.0,1.0);
			double p = dist(generator);
			if (p < prob_stem_cell_divides_to_cancer_cell)
			{
				d.cell_type = CANCER_CELL_TYPE;
			}
		}
		

		return d;
	}

	CellAction Move(vector<Cell>& cell_vec,unsigned char* cell_type_map, const int w, const int h,
			const int num_chem_entities)
	{
		CellAction d(CELL_ACTION_MOVE,0,0,0);
		int ynew = 0;
		int xnew =0;

		pair<int,int> newpos = FindPosToMove(cell_type_map, w, h);


		xnew = newpos.first;
		ynew = newpos.second;

		if ((xnew >0) && (ynew>0))
		{
			cell_type_map[y*w +x] = NULL_CELL_TYPE;
			cell_type_map[ynew*w +xnew] =cell_type;
		}

		d.xpos = xnew;
		d.ypos = ynew;
		d.cell_type = cell_type;
		return d;

	}

	void UpdateActionWeights()
	{

		double avgo2 = 0.0;	
		for (int i=0; i< sensed_chem_vals[O2_INDEX].size(); i++)
		{
			avgo2 += sensed_chem_vals[O2_INDEX][i];
		}
		avgo2 /= sensed_chem_vals[O2_INDEX].size();

		double avg_glucose = 0.0;
		for (int i=0; i< sensed_chem_vals[GLUCOSE_INDEX].size(); i++)
		{
			avg_glucose += sensed_chem_vals[GLUCOSE_INDEX][i];
		}
		avg_glucose /= sensed_chem_vals[GLUCOSE_INDEX].size();

		double avg_vegf = 0.0;
		for (int i=0; i< sensed_chem_vals[VEGF_INDEX].size(); i++)
		{
			avg_vegf += sensed_chem_vals[VEGF_INDEX][i];
		}
		avg_vegf /= sensed_chem_vals[VEGF_INDEX].size();

		double avg_lactic = 0.0;
		for (int i=0; i< sensed_chem_vals[LACTIC_INDEX].size(); i++)
		{
			avg_lactic += sensed_chem_vals[LACTIC_INDEX][i];
		}
		avg_lactic /= sensed_chem_vals[LACTIC_INDEX].size();

		if (cell_type == BLOOD_VESSEL_TYPE)
		{
			//Blood vessels dont die and move. They can only divide
			prob_wt_move =0.0;
			prob_wt_die =0.0;

			//they divide and move depending on the vegef concentration
			//so we only control the epsilon part
			prob_wt_div =1.0;

			//picewise linear
			if (avg_vegf < 10 )
				epsilon = 0.001;
			else if ((avg_vegf >=10.0) && (avg_vegf <= 18.0) )
			{
				epsilon = (18 - avg_vegf)*0.001/(18 -10) + (avg_vegf -10)*0.4/(18-10);
			}
			else
			{
				epsilon = 0.4;
			}
		}
		else if (cell_type == NORMAL_CELL_TYPE)
		{
			//normal cells never move. 
			//die only if very low in glucose
			//dive with low probability 

			//for now let us assume they don't do anything at all

			//TODO - normal to cancer conversion
			prob_wt_move =0.0;
			prob_wt_die =0.0;

			prob_wt_div =0.01;

			epsilon = 1.0;

		}
		else if (cell_type == CANCER_CELL_TYPE)
		{
			//Cancer cells divide, move and die rapidly
			prob_wt_move =0.3;
			prob_wt_die =0.3;
			prob_wt_div =1.7;//0.3;

			//if very low in o2, they can die
			if (avgo2 < 5.0)
			{
				prob_wt_die += 0.1;
			}

			//low on glucose, they can die wih higher prob
			if (avg_glucose < 40)
			{
				prob_wt_die += 0.2;
			}
			//TODO - movement prob depends on local lactic acid concentration also

			epsilon = 0.5;

		}
		else if ( cell_type == CANCER_STEM_CELL_TYPE)
		{
			//Cancer stem cells divide more, move and die rapidly
			prob_wt_move =0.2;
			prob_wt_die =0.2;
			prob_wt_div =1.8;//0.5;

			//if very low in o2, they can die
			if (avgo2 < 5.0)
			{
				prob_wt_die += 0.1;
			}

			//low on glucose, they can die wih higher prob
			if (avg_glucose < 40)
			{
				prob_wt_die += 0.2;
			}

			epsilon = 0.7;

		}
		else
		{
			cout <<"BAD Cell type in Update Weights \n";
		}

	}

	
	CellAction Update(vector<Cell>& cell_vec, 
				unsigned char* cell_type_map, 
				const int num_chem_entities,
				const vector<double*>& conc_array,
				const int w, const int h,
				const int timestep)
	{
		if ( (just_born == true) || (dead_cell == true))
		{
			just_born = false;
			CellAction cell_action(CELL_ACTION_NULL,0,0,cell_type);
			return cell_action;
		}

		//Sense the Environment 
		Sense(conc_array,w, h,timestep);

		//Update Action weights (divide, move, die)
		UpdateActionWeights();

		//Choose an Action and Act (Note: Act will update the cell_vec and cell_type map)
		CellAction cell_action =ChooseActionAndAct(cell_vec, 
												cell_type_map, 
												num_chem_entities,
												conc_array,
												w,h,
												timestep);
		
		return cell_action;

	}

};

class TMESimulator
{
protected:
	int w;
	int h;
	int NSimSteps;
	int NDiffuseStepsPerSimStep;
	int num_chem_entities;
	vector<double*> conc_array; //one w*h sized array for every conc
	vector<double> diff_const;
	unsigned char* cell_type_map; //each location has a marker 0 = no cell, 1= bloodcell, 2 = cancer cell etc
	double* conc_next_temp_array;
	vector<Cell> cell_vec;
	// This will keep changing after every update cell state
	// They are boundaries of various cell types e.g vegf sinks will be just the
	//boundaries of blood vessels. All cells except blood vessels are O2 and Glucose Sinks
	//
	vector<unsigned char*> conc_src_masks;
	vector<unsigned char*> conc_sink_masks;
	vector<double> srcval;
	vector<double> sinkval;
	//Just one chemical entity (single plane of data)
	//assume edges are always a sink
	//c = kd*c +(1-kd)*(sum of neigbhors)
	void DiffuseStep(double* conc_array, double* conc_temp, const int w, const int h,
					 const double kD,//diffusion constant
					 unsigned char* src_mask, unsigned char* sink_mask,unsigned char* cell_type_map, 
					 const double srcval,
					 const double sinkval, const int n_steps =1)

	{
		//cout <<" In Diffuse Step \n";

		double* ca = conc_array;
		double* ct = conc_temp;
		memcpy(ct,ca,w*h*sizeof(double));

		for (int s =0; s< n_steps; s++)
		{
			for (int i=1; i< h-1; i++)
			{
				for (int j =1; j< w-1; j++)
				{
					const double l=ca[i*w+j-1];
					const double r=ca[i*w+j+1];
					const double t=ca[(i-1)*w+j];
					const double b=ca[(i+1)*w+j];

					if (sink_mask[i*w+j])
					{
						ct[i*w+j] = sinkval;
					}
					else if (src_mask[i*w+j])
					{
						ct[i*w+j] = srcval;
					}
					else
					{
						//Diffuse the entity
						
						ct[i*w+j] = kD*(ca[i*w+j]) + ((1-kD)*((l+r+t+b)/4.0));
						// if (ct[i*w+j] >0.0)
						// {
						// 	printf("Diffuse Step: (i,j) = (%d,%d)) c,l,r,t,b =%f,%f,%f,%f,%f\n",i,j,ca[i*w+j],l,r,t,b);
						// }
					}
					//printf("Diffuse Step: (i,j) = (%d,%d)) c,l,r,t,b =%f,%f,%f,%f,%f\n",i,j,ca[i*w+j],l,r,t,b);
					
				}
			}
			
			double* tmp=NULL;
			tmp= ca;
			ca = ct;
			ct = tmp;
		}
		if((n_steps%2)!=0)
		{
			memcpy(ca,ct,w*h*sizeof(double));
		}


	}

	vector< pair<int,int> > create_ellipse(float a, float b, int cx, int cy)
	{
		vector< pair<int, int> > v;
		float tol = 2.0;
		for(int y=cy-b; y<=cy+b; y++) 
		{
    		for(int x=cx-a; x<=cx+a; x++)
    		{
        		if(abs(x*x*a*a+y*y*b*b-a*a*b*b)<=tol)
        		{
        			v.push_back(make_pair(x,y));
        		}
            	
    		}
		}	
		return v;
	}




	void InitCellState(
							const vector< pair<double,double> >& centers,
							const vector<double>& std_dev,
							const int num_cells
							//const long rand_seed = 42
							)
	{
		//find the percentages to mark

		//mark the positions of cancer cells(some as stem cells, some as non-stem cells).
		//mark the position of normal cells.
		//mark the position of immune cells.

		//Generate the normal cells.
		// We assume that the normal cells will be diffused throughout in clumps
		//TODO - ideally they will be arranged around ellipses. Use ellipse code to generate them

		// default_random_engine generator;
		// random_device myRandomDevice;
		// //unsigned long rand_seed = myRandomDevice();
		// unsigned long rand_seed = 42;
		// generator.seed(rand_seed);


		vector<double> cell_prob_vec;
		cell_prob_vec.push_back(prob_cancer_cell);
		cell_prob_vec.push_back(prob_cancer_stem_cell);
		cell_prob_vec.push_back(prob_normal_cell);

		discrete_distribution<size_t> cell_type_distr(cell_prob_vec.begin(), cell_prob_vec.end());

		uniform_int_distribution<> center_dist(0,centers.size()-1);

		for (int d =0; d<num_cells; d++)
		{
			//cout <<"d = " <<d <<"\n";
			const int i = center_dist(generator);
			normal_distribution<double> distribution(0.0,std_dev[i]);
			double vx = distribution(generator);
			double vy = distribution(generator);

			int xpos = int(centers[i].first + vx);
			int ypos = int(centers[i].second +vy);

			if ( (xpos <0) || (xpos > w-2) ||
				 (ypos <0) || (ypos > h-2)
				 )
			{
				continue;
			}

			int cell_type = cell_type_map[ypos*w + xpos];
			if (cell_type == NULL_CELL_TYPE)
			{
				//cell_state_array[ypos*w*nch+ xpos*nch +0] = 1;

				//create the type of cell
				int typ = cell_type_distr(generator);
				unsigned char ctyp = 0;
				if (typ ==0)
				{
					ctyp = CANCER_CELL_TYPE;
				}
				else if ( typ ==1)
				{
					ctyp = CANCER_STEM_CELL_TYPE;
				}
				else if (typ == 2)
				{
					ctyp = NORMAL_CELL_TYPE;
				}
				else
				{
					cout << "ERROR - trying to generate a non existant cell type\n";
				}
				Cell C(xpos, ypos, ctyp,cell_type_map,w,h,num_chem_entities);

				cell_vec.push_back(C);
				cell_type_map[ypos*w + xpos] = ctyp;
			}

		}

		//cout <<"Total Num of Cells = " << cell_vec.size() << "\n";

		//cout <<"Init Blood Vessels ...\n";

		//BLOOD VESSELS
		//Initialize Blood Vessels.
		//We will assume the that the bottom part of the image has one major artery 
		//Sub artieies/capilaries connect the centers to a random point on the major artery by a line
		// Few more random cells are chosen and connected from some existing blood vessel point for realism
		for(int r=1;r<BLOOD_VESSEL_BASE;r++)
		{
			for(int j=1;j<w-1;j++)
			{
				Cell C(j, r, BLOOD_VESSEL_TYPE, cell_type_map,w,h,num_chem_entities);
				cell_vec.push_back(C);
				cell_type_map[r*w + j] = BLOOD_VESSEL_TYPE;
			}
		}

		//making the lines that connect the centres to the blood vessel.
		uniform_int_distribution<> blood_vesel_terminus(1,w-1);
		for(int i=0;i<centers.size();i++)
		{
			const int tx = blood_vesel_terminus(generator);
			const int ty = BLOOD_VESSEL_BASE;

			const int cx = centers[i].first;
			const int cy = centers[i].second;

			// will make sure x0 < x1
			double x0,y0,x1,y1;
			if ( tx < cx)
			{
				x0 = tx; y0 =ty;
				x1 = cx; y1 = cy;
			}
			else 
			{
				x0 = cx; y0 =cy;
				x1 = tx; y1 = ty;
			}

			//draw the line 
			const double m = (y1 -y0)/(x1 -x0);
			const double c = y0 - m*x0;
			for (int sx = x0; sx < x1; sx++)
			{
				const int sy = m*sx+c;
				Cell C(sx, sy, BLOOD_VESSEL_TYPE, cell_type_map,w,h,num_chem_entities);
				cell_vec.push_back(C);
				cell_type_map[sy*w + sx] = BLOOD_VESSEL_TYPE;
			}

		}
		//cout <<"Total Num of Cells + Blood Vessel Cells= " << cell_vec.size() << "\n";


	}

	void UpdateMaps()
	{

		//update src map
		//update sink map

		for (int y =0; y < h; y++)
		{
			for (int x =0; x <w; x++)
			{
				if ((y ==0) || (y ==(h-1)) ||(x ==0) ||(x ==(w-1)))
				{
					//These are permanent sinks for all chem entites
					for (int c = 0; c < num_chem_entities; c++)
					{
						conc_src_masks[c][y*w+x] = 0;
						conc_sink_masks[c][y*w+x] = 1;
						conc_array[c][y*w+x] = sinkval[c];
					}
				}
				else
				{
					const int ctype = cell_type_map[y*w+x];
					if (ctype == NULL_CELL_TYPE)
					{
						for (int c = 0; c < num_chem_entities; c++)
						{
							conc_src_masks[c][y*w+x] = 0;
							conc_sink_masks[c][y*w+x] = 0;
						}
					}
					else 
					{
						//cout <<" found Non Null Cell\n";

						for (int c = 0; c < num_chem_entities; c++)
						{

							if (c == VEGF_INDEX)
							{
								if ( (ctype == CANCER_CELL_TYPE) || 
									 (ctype == CANCER_STEM_CELL_TYPE) 
									 )
								{
									conc_src_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = srcval[c];

									conc_sink_masks[c][y*w+x] = 0;
								}
								else if (ctype == BLOOD_VESSEL_TYPE)
								{
									conc_src_masks[c][y*w+x] = 0;
									
									conc_sink_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = sinkval[c];
								}
								
							}
							else if (c == O2_INDEX)
							{
								if ( (ctype == CANCER_CELL_TYPE) || 
									 (ctype == CANCER_STEM_CELL_TYPE)  ||
									 (ctype == NORMAL_CELL_TYPE)
									 )
								{
									conc_src_masks[c][y*w+x] = 0;

									conc_sink_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = sinkval[c];

								}
								else if (ctype == BLOOD_VESSEL_TYPE)
								{
									conc_src_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = srcval[c];

									conc_sink_masks[c][y*w+x] = 0;
								}
							}
							else if (c == GLUCOSE_INDEX)
							{
								if ( (ctype == CANCER_CELL_TYPE) || 
									 (ctype == CANCER_STEM_CELL_TYPE)  ||
									 (ctype == NORMAL_CELL_TYPE)
									 )
								{
									conc_src_masks[c][y*w+x] = 0;

									conc_sink_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = sinkval[c];
								}
								else if (ctype == BLOOD_VESSEL_TYPE)
								{
									conc_src_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = srcval[c];

									conc_sink_masks[c][y*w+x] = 0;
								}
							}
							else if (c == LACTIC_INDEX)
							{
								if ( (ctype == CANCER_CELL_TYPE) || 
									 (ctype == CANCER_STEM_CELL_TYPE)  
									 )
								{
									conc_src_masks[c][y*w+x] = 1;
									conc_array[c][y*w+x] = srcval[c];

									conc_sink_masks[c][y*w+x] = 0;
								}
							}
							else
							{
								cout <<"BAD value of chemical entity\n";
							}

						}
					}
				}
			}
		}
	}
	

public:
	TMESimulator(int wid, int ht, int ns)
	{
		w=wid+1;
		h=ht+1;
		NSimSteps=ns;
		NDiffuseStepsPerSimStep = 10;//TODO configure this

		num_chem_entities = 4;
		cell_type_map=(unsigned char*)malloc(w*h*sizeof(unsigned char));
		memset(cell_type_map,0,w*h*sizeof(unsigned char));
		conc_next_temp_array=(double*)malloc(w*h*sizeof(double));
		memset(conc_next_temp_array,0,w*h*sizeof(double));


		for(unsigned int k=0; k<num_chem_entities; k++)
		{
			//cout<<"k= "<<k<<"\n";
			double* tca = (double*)malloc(w*h*sizeof(double));
			memset(tca,0,w*h*sizeof(double));
			conc_array.push_back(tca);
			diff_const.push_back(0.1);

			unsigned char* sm1 = (unsigned char*)malloc(w*h*sizeof(unsigned char));
			memset(sm1,0,w*h*sizeof(unsigned char));
			conc_src_masks.push_back(sm1);

			unsigned char* sm2 = (unsigned char*)malloc(w*h*sizeof(unsigned char));
			memset(sm2,0,w*h*sizeof(unsigned char));
			conc_sink_masks.push_back(sm2);

			if (k == O2_INDEX)
			{
				srcval.push_back(BLOOD_VESSEL_O2_CONC);
				sinkval.push_back(0);
			}
			else if ( k == GLUCOSE_INDEX)
			{
				srcval.push_back(BLOOD_VESSEL_GLUCOSE_CONC);
				sinkval.push_back(0);
			}
			else if ( k == VEGF_INDEX)
			{
				srcval.push_back(TUMOR_VEGF_CONC);
				sinkval.push_back(0);
			}
			else if ( k == LACTIC_INDEX)
			{
				srcval.push_back(TUMOR_LACTIC_ACID_CONC);
				sinkval.push_back(0);
			}
			

			//cout <<"num_chem_entities: " << num_chem_entities <<"\n";
			//cout<<" Done alloc for chem entity " <<k <<"\n";
 		}

		const int Numcells = wid*ht*0.3;//TODO Configure

		vector<pair<double,double> > Centres;
		vector<double> Stddevs;
		//Centres.push_back(make_pair(ht/4,wid/4));
		//Centres.push_back(make_pair(ht/3,wid/3));
		Centres.push_back(make_pair(ht/2,wid/2));
		Stddevs.push_back(4);
		//Stddevs.push_back(3);
		//Stddevs.push_back(2);

		//cout << "Starting Update Cell State .... \n";

		InitCellState(Centres, Stddevs, Numcells);
		//cout << "Update Cell State Done \n";
		UpdateMaps();
		//cout << "Update Maps Done \n";


	}

	void UpdateConc(const int ns =1)
	{
		for ( int s =0; s <ns; s++)
		{
			for (int i=0; i < conc_array.size(); i++)
			{
				DiffuseStep(conc_array[i],conc_next_temp_array,w,h,
					diff_const[i],
					conc_src_masks[i],
					conc_sink_masks[i], cell_type_map, srcval[i] , sinkval[i],NDiffuseStepsPerSimStep);

				//cout<<"Diffuse step "<<s<< " Conc id " <<i <<" done\n";

			}
			
		}
	};


	void SaveImg(double* diffarray,int w, int h,char* filename)
	{

	}

	double AvgConc(const double* conc_array, const int w, const int h)
	{
		double avg = 0.0;
		for (int y=1; y< h-1; y++)
		{
			for (int x =0; x< w-1; x++)
			{
				//cout <<"("<<x<<","<<y<<")"<<" " << conc_array[y*w +x] <<"\n";
				avg += conc_array[y*w +x];
			}
		}
		//cout <<"\n";
		//cout <<"AvgConc - w = " << w << " h = " <<h <<" avg = " << avg <<"\n";
		return avg/((w-1)*(h-1));
	}

	void Measure(const int timestep)
	{
		//bounding box of cancer cells
		//bbox of normal cells
		//number of cancer cells, number of cancer stem cells
		//number of blood vessel cells
		//avg conc of all four chemicals
		int csx0=w,csy0=h,csx1=-1,csy1=-1;
		int cx0=w,cy0=h,cx1=-1,cy1=-1;
		int nx0=w,ny0=h,nx1=-1,ny1=-1;
		int bx0=w,by0=h,bx1=-1,by1=-1;
		int num_cancer_cells=0, num_normal_cells=0,num_blood_cells=0,num_cancer_stem_cells=0;
		int num_dead_cells=0,num_movements=0;
		for (int i =0; i<cell_vec.size(); i++)
		{
			if (cell_vec[i].dead_cell == false)
			{
				if (cell_vec[i].cell_type == NORMAL_CELL_TYPE)
				{
					if (cell_vec[i].x < nx0)
					{
						nx0 = cell_vec[i].x;
					}
					if (cell_vec[i].y < ny0)
					{
						ny0 = cell_vec[i].y;
					}
					if (cell_vec[i].x > nx1)
					{
						nx1 = cell_vec[i].x;
					}
					if (cell_vec[i].y > ny1)
					{
						ny1 = cell_vec[i].y;
					}
					num_normal_cells += 1;
				}
				else if (cell_vec[i].cell_type == CANCER_CELL_TYPE)
				{
					if (cell_vec[i].x < cx0)
					{
						cx0 = cell_vec[i].x;
					}
					if (cell_vec[i].y < cy0)
					{
						cy0 = cell_vec[i].y;
					}
					if (cell_vec[i].x > cx1)
					{
						cx1 = cell_vec[i].x;
					}
					if (cell_vec[i].y > cy1)
					{
						cy1 = cell_vec[i].y;
					}
					num_cancer_cells += 1;
				}
				else if (cell_vec[i].cell_type == CANCER_STEM_CELL_TYPE)
				{
					if (cell_vec[i].x < csx0)
					{
						csx0 = cell_vec[i].x;
					}
					if (cell_vec[i].y < csy0)
					{
						csy0 = cell_vec[i].y;
					}
					if (cell_vec[i].x > csx1)
					{
						csx1 = cell_vec[i].x;
					}
					if (cell_vec[i].y > csy1)
					{
						csy1 = cell_vec[i].y;
					}
					num_cancer_stem_cells += 1;
				}
				else if (cell_vec[i].cell_type == BLOOD_VESSEL_TYPE)
				{
					if (cell_vec[i].x < bx0)
					{
						bx0 = cell_vec[i].x;
					}
					if (cell_vec[i].y < by0)
					{
						by0 = cell_vec[i].y;
					}
					if (cell_vec[i].x > bx1)
					{
						bx1 = cell_vec[i].x;
					}
					if (cell_vec[i].y > by1)
					{
						by1 = cell_vec[i].y;
					}
					num_blood_cells += 1;
				}
			}
			else
			{
				num_dead_cells += 1;
			}
		}
		double avg_glucose = AvgConc(conc_array[GLUCOSE_INDEX], w, h);
		double avg_o2 = AvgConc(conc_array[O2_INDEX], w, h);
		double avg_lactic = AvgConc(conc_array[LACTIC_INDEX], w, h);
		double avg_vegf = AvgConc(conc_array[VEGF_INDEX], w, h);

		double area_cancer_stem_cells = 0.0;
		if (( csx1 >= csx0) && (csy1 >= csy0))
			area_cancer_stem_cells = (csy1+1 - csy0)*(csx1+1 - csx0);

		double area_cancer_cells = 0.0;
		if (( cx1 >= cx0) && (cy1 >= cy0))
			area_cancer_cells = (cy1+1 - cy0)*(cx1+1 - cx0);

		double area_normal_cells = 0.0;
		if (( nx1 >= nx0) && (ny1 >= ny0))
			area_normal_cells = (ny1+1 - ny0)*(nx1+1 - nx0);

		double area_blood_cells = 0.0;
		if (( bx1 >= bx0) && (by1 >= by0))
			area_blood_cells = (by1+1 - by0)*(bx1+1 - bx0);

		cout << timestep <<",";
		cout << csx0 <<"," << csy0 <<"," << csx1 <<"," << csy1 <<", ," ;
		cout << cx0 <<"," << cy0 <<"," << cx1 <<"," << cy1 <<", ,";
		cout << nx0 <<"," << ny0 <<"," << nx1 <<"," << ny1 <<", ,";
		cout << bx0 <<"," << by0 <<"," << bx1 <<"," << by1 <<", ,";
		cout << area_cancer_stem_cells <<"," << area_cancer_cells <<"," << area_normal_cells <<"," << area_blood_cells <<", ,";
		cout << num_cancer_stem_cells <<"," << num_cancer_cells << "," << num_normal_cells <<"," << num_blood_cells<<", ,";
		cout <<avg_vegf <<"," << avg_lactic <<"," << avg_o2 << "," << avg_glucose;
		cout <<"\n";
		
		//TODO add number of movements also
		return;

	}


	void UpdateCellState(const int ts)
	{
		vector<int> dead_cell_indices;
		vector <int> moved_cell_indices;
		vector<int> no_action_cell_indices;
		vector <CellAction> new_cells;
		for (int i =0; i<cell_vec.size(); i++)
		{
			//each cell will do something depending on it's type
			CellAction cell_action = cell_vec[i].Update(cell_vec, 
														cell_type_map, 
														num_chem_entities,
														conc_array,
														w,h,ts);
			if (cell_action.action_type == CELL_ACTION_DIE)
			{
				//mark this for removal
				dead_cell_indices.push_back(i);
				cell_type_map[cell_vec[i].y*w +cell_vec[i].x] = NULL_CELL_TYPE;
			}
			else if (cell_action.action_type == CELL_ACTION_MOVE)
			{
				if ( (cell_action.xpos > 0) && (cell_action.ypos >0))
				{
					moved_cell_indices.push_back(i);
					cell_type_map[cell_vec[i].y*w +cell_vec[i].x] = NULL_CELL_TYPE;
					cell_type_map[cell_action.ypos*w +cell_action.xpos] = cell_action.cell_type;
				}

			}
			else if (cell_action.action_type == CELL_ACTION_DIVIDE)
			{
				if ( (cell_action.xpos > 0) && (cell_action.ypos >0))
				{
					new_cells.push_back(cell_action);
					cell_type_map[cell_action.ypos*w +cell_action.xpos] = cell_action.cell_type;
				}
			}
			else if (cell_action.action_type == CELL_ACTION_NULL)
			{
				no_action_cell_indices.push_back(i);
			}
			else
			{
				cout <<" BAD CELL action\n";
			}
		}

		// Add all the new cells
		for (int i =0; i< new_cells.size(); i++)
		{
			CellAction cell_action = new_cells[i];
			cell_vec.push_back(Cell(cell_action.xpos, cell_action.ypos, cell_action.cell_type,
										 cell_type_map, w,  h, num_chem_entities));
			cell_type_map[cell_action.ypos*w +cell_action.xpos] = cell_action.cell_type;
		}

		UpdateMaps();

		//cout <<"nCells Added = " << new_cells.size() <<" Dead = " << dead_cell_indices.size() << "\n";
	}

};

int main(int argc, char* argv[])
{
	//cout <<"usage: simsim imageDim nSimSteps";

	//hello<<<1,1>>>();
	const int N = atoi(argv[1]);
	//cout <<"N =" <<N <<"\n";

	const int NSimSteps = atoi(argv[2]);

	const int n_chem_entities = 4; //VEGF, O2, Ph, Glucose, 

	//Initialize the global random number generator
	random_device myRandomDevice;
	//unsigned long rand_seed = myRandomDevice();
	unsigned long rand_seed = 42;
	generator.seed(rand_seed);

	TMESimulator T(N, N, NSimSteps);

	cout << "timestep" <<",";
	cout << "csx0" <<"," << "csy0" <<"," << "csx1" <<"," << "csy1" <<", ,";
	cout << "cx0" <<"," << "cy0" <<"," << "cx1" <<"," << "cy1" <<", ,";
	cout << "nx0" <<"," << "ny0" <<"," << "nx1" <<"," << "ny1" <<", ,";
	cout << "bx0" <<"," << "by0" <<"," << "bx1" <<"," << "by1" <<", ,";
	cout << "area_cancer_stem_cells" <<"," << "area_cancer_cells" <<"," << "area_normal_cells" <<"," << "area_blood_cells" <<", ,";
	cout << "num_cancer_stem_cells" <<"," << "num_cancer_cells"<< "," << "num_normal_cells" <<"," << "num_blood_cells"<<", ,";
	
	cout <<"avg_vegf" <<"," << "avg_lactic" <<"," << "avg_o2" << "," << "avg_glucose";
	cout << "\n";

	T.Measure(0);

	T.UpdateConc(1000);
	T.Measure(0);

	for (int i =0; i <NSimSteps; i++)
	{
		//cout <<"Nsimstep = " << i << "\n";
		T.UpdateCellState(i);
		T.UpdateConc(10);
		//Measure the obsevables and store 
		T.Measure(i+1);
	}
	//cout <<" All Sim Steps Done!!!\n";
	//Draw the graphs

	return 0;
	
}