// ************************************************************************************* //
//    Program to solve the Facility Location Reformulation of the General Capacitated    //
//             Lot-Sizing Problem with Multiple Storage Locations (RCLSP-MSL)            //
//                         by a Sequential Heuristic - Version 3                         //
//	 																					 //
//	  Used in Gislaine Mara Melega Pos-doctoral											 //
//    Copyright 2020 - date:  04/2022													 //
// ************************************************************************************* //



//Libraries
//#include <stdafx.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloexpression.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <cstdlib>
#include <vector>



//macro necessary for portability
ILOSTLBEGIN	



// ********************************************************************************************** //
// ******************************** BEGINING OF THE MAIN PROGRAM ******************************** //
// ********************************************************************************************** //


int main(int argc, char **argv)
{
 
    //Input Data File
    ifstream in(argv[1]);


    //Output Data File 1
    ofstream out(argv[2]);


    //Problem enviroment: env
    IloEnv env;


    try{  


          int  t, i, l,													  //indexes to: time periods, items, locations
			 
			   j, k, tau;												  //other indexes

		  
		  //Variables name
		  char setupitem[15];											  //setup of item
		  char proditem[15];											  //production of item
		  char stockitemloc[15];										  //stock of item in location

		  char setuplocation[15];									      //setuplocation
		  char inflowitem[15];											  //inflow of item in location 
		  char outflowitem[15];											  //outflow of item in location
		  char assignitemlocation[15];									  //assignment of item to storage location
		  char relocationitem[15];										  //relocation of item between locations

		  char facilityref[15];											  //facility location reformulation



// ************************************ // 
//    Creation of the Parameters and    //
//     Data Reading from Input File     // 
// ************************************ //


		  //Indexes
		  //Time periods
		  IloInt T;
		  
		  //Items
		  IloInt I;

		  //Locations
		  IloInt L;



		  //Read only the indexes
		  //to create the parameters
		  if (in) {

				     in  >>  T;
				     in  >>  I;
				     in  >>  L;

		  }
	        else {
                    cerr << "No such file: " << "dataRCLSPMSL.dat" << endl;
                    throw(1);
			}

          // *************************************************************


		  //Production of Items
          //Setup cost of item 
		  IloNumArray sc(env, I);										

		  //Production cost per unit of item
          IloNumArray vc(env, I);

		  //Inventory holding cost per unit of item
		  IloNumArray hc(env, I);

		  //Production consumption of capacity per unit of item 
		  IloNumArray vt(env, I);

		  //Production capacity in each period
		  IloNumArray Cap(env, T);

		  //Demand to be fulfilled of item i in period t
		  IloArray<IloNumArray> d(env, I);
		  for(i=0; i<I; i++)
		     d[i] = IloNumArray(env, T);					


		  //Inventory of items
		  //Fixed cost of having a positive
		  //inventory level at location l
		  IloNumArray g(env, L);

		  //Unit handling cost of item i at location l
		  IloArray<IloNumArray> ha(env, I);
		  for(i=0; i<I; i++)
		     ha[i] = IloNumArray(env, L);	

		  //Consumption of storage capacity per unit of item
		  IloNumArray cs(env, I);

		  //Storage capacity of each location in each period
		  IloNumArray H(env, L);

		  //Compatibility between item i and locations l
		  IloArray<IloNumArray> alpha(env, I);
		  for(i=0; i<I; i++)
		     alpha[i] = IloNumArray(env, L);	

		  //Compatibility between item i an j
		  IloArray<IloNumArray> beta(env, I);
		  for(i=0; i<I; i++)
		     beta[i] = IloNumArray(env, I);	

		  //Unit moving cost of item i
		  //from location k to location l 
		  IloArray<IloArray<IloNumArray> > r(env, I);
		  for(i=0; i<I; i++){
		     r[i] = IloArray<IloNumArray> (env, L); 
		     for(l=0; l<L; l++){
		        r[i][l] = IloNumArray(env, L);
			 }								
		  }          


		  //Big M
		  IloInt BigM;



		  //Read the remaining parameters
		  if (in) {
			          
				     for(t=0; t<T; t++)
					    in >> Cap[t];

					 
					 for(i=0; i<I; i++){
 					    in >> vc[i];
						in >> sc[i];
						in >> hc[i];
						in >> vt[i];
						in >> cs[i];
					 }


					 for(i=0; i<I; i++)
					    for(t=0; t<T; t++)
						   in >> d[i][t];


					 for(l=0; l<L; l++){
					    in >> H[l];
						in >> g[l];
					 }


					 for(i=0; i<I; i++)
				        for(l=0; l<L; l++)
						   in >> ha[i][l];


					 for(i=0; i<I; i++)
				        for(l=0; l<L; l++)
						   in >> alpha[i][l];


					 for(i=0; i<I; i++)
                        for(j=0; j<I; j++)
						   in >> beta[i][j];				   


					 for(i=0; i<I; i++)
					    for(l=0; l<L; l++)
						   for(k=0; k<L; k++)
						      in >> r[i][l][k];


          }//end if(in)
	        else {
                    cerr << "No such file: " << "dataRCLSPMSL.dat" << endl;
                    throw(1);
			}

// ***********************************************************************


	 


// ************************ // 
//    Create the Problem    // 
// ************************ //


		  //Problem
		  IloModel Pmodel(env);


		  //All the variables are considered linear (relaxed values)

		  //Setup of item i in period t
		  IloArray<IloNumVarArray> Y(env, I);
		  for(i=0; i<I; i++){
		     Y[i] = IloNumVarArray(env, T);
			 for(t=0; t<T; t++){
                sprintf(setupitem, "Y_%d_%d", i, t);
				Y[i][t] = IloNumVar(env, 0, 1, setupitem);
			 }
		  }



		  //Production of item i in period t
		  IloArray<IloNumVarArray> X(env, I);
		  for(i=0; i<I; i++){
		     X[i] = IloNumVarArray(env, T);
			 for(t=0; t<T; t++){
                sprintf(proditem, "X_%d_%d", i, t);
				X[i][t] = IloNumVar(env, 0, IloInfinity, proditem);
			 }
		  }



		  //Inventory of item i at location l in period t
		  IloArray<IloArray<IloNumVarArray> > S(env, I);
		  for(i=0; i<I; i++){
		     S[i] = IloArray<IloNumVarArray> (env, L); 
		     for(l=0; l<L; l++){ 
		        S[i][l] = IloNumVarArray(env, T);								
				for(t=0; t<T; t++){ 
				   sprintf(stockitemloc, "S_%d_%d_%d", i, l, t);
                   S[i][l][t] = IloNumVar(env, 0, IloInfinity, stockitemloc);
				}
			 }
		  }



		  //Use of location l in period t
		  //i.e., there is a positive inventory
		  //in location l at the end of period t
		  IloArray<IloNumVarArray> Z(env, L);
		  for(l=0; l<L; l++){
		     Z[l] = IloNumVarArray(env, T);
			 for(t=0; t<T; t++){
                sprintf(setuplocation, "Z_%d_%d", l, t);
				Z[l][t] = IloNumVar(env, 0, 1, setuplocation);
			 }
		  }



		  //Inflow of item i at location l in period t
		  IloArray<IloArray<IloNumVarArray> > Dp(env, I);
		  for(i=0; i<I; i++){
		     Dp[i] = IloArray<IloNumVarArray> (env, L); 
		     for(l=0; l<L; l++){ 
		        Dp[i][l] = IloNumVarArray(env, T);								
				for(t=0; t<T; t++){ 
				   sprintf(inflowitem, "Dp_%d_%d_%d", i, l, t);
                   Dp[i][l][t] = IloNumVar(env, 0, IloInfinity, inflowitem);
				}
			 }
		  }



		  //Outflow of item i at location l in period t
		  IloArray<IloArray<IloNumVarArray> > Dm(env, I);
		  for(i=0; i<I; i++){
		     Dm[i] = IloArray<IloNumVarArray> (env, L); 
		     for(l=0; l<L; l++){ 
		        Dm[i][l] = IloNumVarArray(env, T);								
				for(t=0; t<T; t++){ 
				   sprintf(outflowitem, "Dm_%d_%d_%d", i, l, t);
                   Dm[i][l][t] = IloNumVar(env, 0, IloInfinity, outflowitem);
				}
			 }
		  }



		  //Assignment of item i at storage location l in period t
		  IloArray<IloArray<IloNumVarArray> > W(env, I);
		  for(i=0; i<I; i++){
		     W[i] = IloArray<IloNumVarArray> (env, L); 
		     for(l=0; l<L; l++){ 
		        W[i][l] = IloNumVarArray(env, T);								
				for(t=0; t<T; t++){ 
				   sprintf(assignitemlocation, "W_%d_%d_%d", i, l, t);
                   W[i][l][t] = IloNumVar(env, 0, 1, assignitemlocation);
				}
			 }
		  }



		  //Relocation of item i from location l to location k in period t
		  IloArray<IloArray<IloArray<IloNumVarArray> > > V(env, I);
		  for(i=0; i<I; i++){
		     V[i] = IloArray<IloArray<IloNumVarArray> > (env, L); 
			 for(l=0; l<L; l++){
				 V[i][l] = IloArray<IloNumVarArray> (env, L); 
				 for(k=0; k<L; k++){ 
					V[i][l][k] = IloNumVarArray (env, T);
					for(t=0; t<T; t++){
					   sprintf(relocationitem, "V_%d_%d_%d_%d", i, l, k, t);
					   V[i][l][k][t] = IloNumVar(env, 0, IloInfinity, relocationitem);
					}
				 }
			  }
		  }



		  //Facility location reformulation 
		  //number of units of item i produced in period t to
		  //meet the demand in a posterior period tau, tau >= t
		  IloArray<IloArray<IloNumVarArray> > FL(env, I);
		  for(i=0; i<I; i++){
		     FL[i] = IloArray<IloNumVarArray> (env, T); 
		     for(t=0; t<T; t++){ 
		        FL[i][t] = IloNumVarArray(env, T);								
				for(tau = t; tau<T; tau++){ 
				   sprintf(facilityref, "FL_%d_%d_%d", i, t, tau);
                   FL[i][t][tau] = IloNumVar(env, 0, IloInfinity, facilityref);
				}
			 }
		  }

		  // *************************************************************



		  //Objective Function
		  IloExpr objective(env);


		  //Costs of setup of items
		  for(t=0; t<T; t++)
			 for(i=0; i<I; i++)
		        objective += sc[i]*Y[i][t];


		  //Costs of production
		  for(t=0; t<T; t++)
			 for(tau = t; tau<T; tau++)
			    for(i=0; i<I; i++)
		           objective += vc[i]*FL[i][t][tau];


		  //Cost of inventory and handling of items at storage locations
		  for(t=0; t<T; t++)
			 for(l=0; l<L; l++)
			    for(i=0; i<I; i++)
		           objective += (hc[i]*S[i][l][t] + ha[i][l]*Dp[i][l][t]);


		  //Costs of using locations
		  for(t=0; t<T; t++)
		     for(l=0; l<L; l++)
			    objective += g[l]*Z[l][t];


		  //Cost of relocation of items between locations
		  for(t=0; t<T; t++)
			 for(l=0; l<L; l++)
				for(k=0; k<L; k++)
			       for(i=0; i<I; i++)
			          objective += r[i][l][k]*V[i][l][k][t];


		  //Problem objective function environment
		  IloObjective Pof = IloMinimize(env, objective);


          //Add the objective function (Pof) to the problem
		  Pmodel.add(Pof);
          objective.end();		     //delet the expression

		  // *************************************************************



		  //InflowOutflow1 constraints environment
		  IloArray<IloRangeArray> InflowOutflow1(env, I);
		  for(i=0; i<I; i++)
		     InflowOutflow1[i] = IloRangeArray(env, T);

		  //InflowOutflow1 constraints
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++){
				IloExpr flow(env);
				
				for(tau=0; tau <= t; tau++)
				   flow += FL[i][tau][t];

			    InflowOutflow1[i][t] = (flow == d[i][t]);
				flow.end();
			 }

		  //Add the InflowOutflow1 constraints to the problem
		  for(i=0; i<I; i++){
             InflowOutflow1[i].setNames("InflowOutflow1");
	         Pmodel.add(InflowOutflow1[i]);
		  }





		  //InflowOutflow2 constraints environment
		  IloArray<IloRangeArray> InflowOutflow2(env, I);
		  for(i=0; i<I; i++)
		     InflowOutflow2[i] = IloRangeArray(env, T);

		  //InflowOutflow2 constraints
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++){
				IloExpr flow(env);
				
				for(tau = t; tau<T; tau++)
				   flow += FL[i][t][tau];

				flow -= d[i][t];

				for(l=0; l<L; l++)
				   flow -= Dp[i][l][t] - Dm[i][l][t];

			    InflowOutflow2[i][t] = (flow == 0);
				flow.end();
			 }

		  //Add the InflowOutflow2 constraints to the problem
		  for(i=0; i<I; i++){
             InflowOutflow2[i].setNames("InflowOutflow2");
	         Pmodel.add(InflowOutflow2[i]);
		  }




		  
		  //BalanceLocation constraint environment
		  IloArray<IloArray<IloRangeArray> > BalanceLocation(env, I);
		  for(i=0; i<I; i++){
		     BalanceLocation[i] = IloArray<IloRangeArray>(env,L);
		     for(l=0; l<L; l++){
		        BalanceLocation[i][l] = IloRangeArray(env, T);
			 }
		  }

		  //BalanceLocation constraint
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++)
		        for(t=0; t<T; t++){
		           IloExpr balance(env);

				   if (t == 0) {
								  balance += Dp[i][l][t] - S[i][l][t] - Dm[i][l][t];

								  for(k=0; k<L; k++)
									 balance += V[i][k][l][t] - V[i][l][k][t];

								  BalanceLocation[i][l][t] = (balance == 0);
								  balance.end();
					}
					  else {

								  balance += S[i][l][t-1] + Dp[i][l][t] - S[i][l][t] - Dm[i][l][t];

								  for(k=0; k<L; k++)
								     balance += V[i][k][l][t] - V[i][l][k][t];

								  BalanceLocation[i][l][t] = (balance == 0);
								  balance.end();
					  }
				}

		  //Add the BalanceLocation constraint to the problem
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++){
                BalanceLocation[i][l].setNames("BalanceLocation");
	            Pmodel.add(BalanceLocation[i][l]);
			 }





		  //Setup constraints environment
		  IloArray<IloArray<IloRangeArray> > Setup(env, I);
		  for(i=0; i<I; i++){
		     Setup[i] = IloArray<IloRangeArray>(env, T);
		     for(t=0; t<T; t++){
		        Setup[i][t] = IloRangeArray(env, T);
			 }
		  }

		  //Setup constraints
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++)
			    for(tau = t; tau<T; tau++){
			       IloExpr set(env);
				
                   set += FL[i][t][tau] - d[i][tau]*Y[i][t];

			    Setup[i][t][tau] = (set <= 0);
				set.end();
		     }

		  //Add the Setup constraints to the problem
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++)
			    for(tau = t; tau<T; tau++){
                   Setup[i][t][tau].setName("Setup");
	               Pmodel.add(Setup[i][t][tau]);
				}


		  


		  //Capacity constraints environment
		  IloRangeArray Capacity(env, T);

		  //Capacity constraints
		  for(t=0; t<T; t++){
		     IloExpr cap(env);

			 for(i=0; i<I; i++)
			    for(tau = t; tau<T; tau++)
			       cap += vt[i]*FL[i][t][tau]; 

			 Capacity[t] = (cap <= Cap[t]);
			 cap.end();
		  }

		  //Add the Capacity constraints to the problem
          Capacity.setNames("Capacity");
	      Pmodel.add(Capacity);





		  //InvAlloc constraints enviroment
		  IloArray<IloArray<IloRangeArray> > InvAlloc(env, I);
		  for(i=0; i<I; i++){
		     InvAlloc[i] = IloArray<IloRangeArray>(env,L);
		     for(l=0; l<L; l++){
		        InvAlloc[i][l] = IloRangeArray(env, T);
			 }
		  }

		  //InvAlloc constraints
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++)
				for(t=0; t<T; t++){
				   IloExpr alloc(env);

				   //Calculate the bigM 
				   int sumd = 0; 
				   double remaincapL = 0;

				   //Sum of the demand
				   for(tau=t; tau<T; tau++)
				      sumd += d[i][tau];

				   //Remaning storage capacity 
				   remaincapL = (H[l]/cs[i]);

				   //Choose the smallest value between the sum
				   //of the demand and remaning storage capacity 
				   if (sumd <= remaincapL) BigM = sumd;
				     else BigM = remaincapL;


                   alloc += S[i][l][t] - BigM*W[i][l][t];

			       InvAlloc[i][l][t] = (alloc <= 0);
				   alloc.end();
				}

		  //Add the InvAlloc constraints to the problem
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++){
                InvAlloc[i][l].setNames("InvAlloc");
	            Pmodel.add(InvAlloc[i][l]);
			 }





		  //CapacityStorage constraints environment
    	  IloArray<IloRangeArray> CapacityStorage(env, L); 
	      for(l=0; l<L; l++){
		     CapacityStorage[l] = IloRangeArray(env, T);
		  }

		  //CapacityStorage
		  for(l=0; l<L; l++)
		     for(t=0; t<T; t++){
			    IloExpr cap(env);

				for(i=0; i<I; i++)
				   cap += cs[i]*S[i][l][t];

				cap -= H[l]*Z[l][t];
			    
				CapacityStorage[l][t] = (cap <= 0);
				cap.end();
			 }

		  //Add the CapacityStorage constraints to the problem
		  for(l=0; l<L; l++){
             CapacityStorage[l].setNames("CapacityStorage");
	         Pmodel.add(CapacityStorage[l]);
		  }





		  //ItemLocatCompat constraints enviroment
		  IloArray<IloArray<IloRangeArray> > ItemLocatCompat(env, I);
		  for(i=0; i<I; i++){
		     ItemLocatCompat[i] = IloArray<IloRangeArray>(env,L);
		     for(l=0; l<L; l++){
		        ItemLocatCompat[i][l] = IloRangeArray(env, T);
			 }
		  }

		  //ItemLocatCompat constraints
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++)
				for(t=0; t<T; t++){
				   IloExpr comp(env);

                   comp += W[i][l][t];

			       ItemLocatCompat[i][l][t] = (comp <= alpha[i][l]);
				   comp.end();
			 }

		  //Add the ItemLocatCompat constraints to the problem
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++){
                ItemLocatCompat[i][l].setNames("ItemLocatCompat");
	            Pmodel.add(ItemLocatCompat[i][l]);
			 }





		  //ItemItemCompat constraints enviroment
		  IloArray<IloArray<IloArray<IloRangeArray> > > ItemItemCompat(env, I);
		  for(i=0; i<I; i++){
             ItemItemCompat[i] = IloArray<IloArray<IloRangeArray> >(env, I);
             for(j=i; j<I; j++){
		        ItemItemCompat[i][j] = IloArray<IloRangeArray>(env,L);
		        for(l=0; l<L; l++){
		           ItemItemCompat[i][j][l] = IloRangeArray(env, T);
				}
			 }
		  }

		  //ItemItemCompat constraints
		  for(i=0; i<I; i++)
			 for(j=i; j<I; j++)
		        for(l=0; l<L; l++)
				   for(t=0; t<T; t++){
				      IloExpr comp(env);

                      comp += W[i][l][t] + W[j][l][t];

			          ItemItemCompat[i][j][l][t] = (comp <= beta[i][j] + 1);
				      comp.end();
				   }

		  //Add the ItemItemCompat constraints to the problem
		  for(i=0; i<I; i++)
			 for(j=i; j<I; j++)
		        for(l=0; l<L; l++){
                   ItemItemCompat[i][j][l].setNames("ItemItemCompat");
	               Pmodel.add(ItemItemCompat[i][j][l]);
				}

		  // *************************************************************


		  //Define CPLEX environment to the problem
          IloCplex Pcplex(Pmodel);


// ***********************************************************************





// ************************************************************************* //
//    Solve the Linear Relaxation of the Reformulated General Capacitated    //
//       Lot-Sizing Problem with Multiple Storage Locations (RCLSP-MSL)      //
//                         by an Optimization Package                        //
// ************************************************************************* //


		  //Objective function value and computational
		  //time of the linear relaxation problem
		  double  OF_LR_RCLSPMSL = 0, Time_LR_RCLSPMSL = 0;



          // ***** Solve the LR_RCLSP-MSL ********************************

		  //Add CPLEX Options
          //Print the output and warnings 
		  //of cplex in the output file
          Pcplex.setOut(out);
          Pcplex.setWarning(out);


		  //Limite the number of threads in
		  //the solution of the linear relation problem
		  Pcplex.setParam(IloCplex::Threads, 1);



		  //Extract the model
		  //Pcplex.extract(Pmodel);
		  //Export the model
		  //Pcplex.exportModel("model.lp" );

	 
		  
          // ****************************************************************************************************
		  out << "***** The Facility Location Reformulation of the General Capacitated Lot-Sizing *****" << endl;
		  out << "*****            Problem with Multiple Storage Locations (RCLSP-MSL)            *****" << endl;       
		  out << "*****                      solved by a Sequential Heuristic                     *****" << endl;
		  out << "*****                                 Version 3                                 *****" << endl;
		  out << endl << endl;
		  out << "* Lower Bound: Linear Relaxation" << endl;
		  out << "* Upper Bound: Sequential Heuristic - Version 3" << endl;
		  out << endl << endl << endl << endl << endl;
		  out << "*********** Lower Bound: Solve the LR_RCLSP-MSL by an Optimization Package **********" << endl;
		  out << endl << endl;


		  //SOLVE the LR_RCLSP-MSL
          Pcplex.solve();


		  out << endl << endl << endl;
		  // ****************************************************************************************************	  





		  // ***** Check the Solution Status *****************************

		  
		  //Check if the solution of the Linear Relaxation Problem is
		  //Optimal-1 or Feasible-2
		  if ((Pcplex.getStatus() == IloAlgorithm::Optimal) || 
			  (Pcplex.getStatus() == IloAlgorithm::Feasible)) {
	

					//Recover the objective function value
				    //of the linear relation problem
					OF_LR_RCLSPMSL = Pcplex.getValue(Pof);


					//Recover the time to solve 
					//the linear relaxation problem
					Time_LR_RCLSPMSL = Pcplex.getTime();

					

					// ****************************************************************************************************
					//Print in the output file
					out << "***** Solution to the Linear Relaxation of the Reformulated General Capacitated *****" << endl;
					out << "*****             Lot-Sizing Problem with Multiple Storage Locations            *****" << endl;
					out << endl << endl;
					out << "Solution Status LR = " << Pcplex.getStatus() << endl;
					out << "Objective Function Value LR = " << OF_LR_RCLSPMSL << endl; 
					out << "Time LR = " << Time_LR_RCLSPMSL << endl;	 
					out << endl << endl;
					out << "*************************************************************************************" << endl;
					out << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl;
					// ****************************************************************************************************


		  }//end if Optimal-1/Feasible-2

		    else {

				    // ****************************************************************************************************
					//Print in the output file
					out << "***** Solution to the Linear Relaxation of the Reformulated General Capacitated *****" << endl;
					out << "*****             Lot-Sizing Problem with Multiple Storage Locations            *****" << endl;
					out << endl << endl;
					out << "Solution Status LR = " << Pcplex.getStatus() << endl;
					out << "NO Solution to the Linear Relaxation Problem" << endl;
					out << endl << endl;
					out << "*************************************************************************************" << endl;
					out << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl;
					// ****************************************************************************************************


			}//end else


// ***********************************************************************





// ******************************************************************* //
//    Solve the Reformulated General Capacitated Lot-Sizing Problem    // 
//			   with Multiple Storage Locations (RCLSP-MSL)			   //
//                by a Sequential Heuristic - Version 3                //
// ******************************************************************* //


//* Sequential Heuristic Problem 1 (SHP1): facility location reformulation of the
//										   capacitated lot-sizing problem (RCLSP)
//* Sequential Heuristic Problem 2 (SHP2): facility location reformulation of the
//										   capacitated lot-sizing problem with multiple storage locations (RCLSP-MSL)

//Note: In SHP2, the CLSP-MSL results in the multiple storage location 
//problem (MSLP), after the fixing of the decision variables from SHP1


		  //Objective function value,
		  //gap, and computational time
		  //of the Sequential Heuristic
		  double OF_SH = IloInfinity, Gap_SH = 0, Time_SH = 0;


		  //Other calculation of a solution
		  //of the Sequential Heuristic
		  double CSetupItem = 0, CInventItem = 0, CHandItem = 0, CSetupLocal = 0, CRelocItem = 0,
				 NSetupItem = 0, NInventItem = 0, NHandItem = 0, NLocalUsed = 0, NRelocItem = 0,
				 TotalOpenSpace = 0, TotalUsedSpace = 0, PercUsedSpace = 0;



		  //Time limit to the
		  //Sequential Heuristic			  
		  double TimeLimit_SH = 1800;

		  
		  //Iterations of the
		  //Sequencial Heuristic 
		  int it_SH = 0;


		  //Objective function value between two
		  //consecutive iterations of the Sequencial Heuristic 
		  double OF_SH_it1 = IloInfinity, OF_SH_it2 = IloInfinity;





		  //Estimation for the handling cost
		  IloNum InitialEpsilon = 0.25;
		  IloNum Epsilon = 0.25;


		  //Create the parameter to the 
		  //estimation for the handling cost
		  IloArray<IloNumArray> est_ha(env, I);
		  for(i=0; i<I; i++)
		     est_ha[i] = IloNumArray(env, L);

		  //First iteration
		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++)
		        est_ha[i][l] = InitialEpsilon*ha[i][l];





		  // ****************************************************************************************************
		  //Print in the output file
		  out << "*********** Upper Bound: Beginning of the Sequential Heuristic - Version 3 **********" << endl;
		  out << endl << endl;
		  out << "* Sequential Heuristic: SHP1 - facility location reformulation of the" << endl;
		  out << "        Version 3              capacitated lot-sizing problem (RCLSP)" << endl;
          out << "                        SHP2 - facility location reformulation of the" << endl;
		  out << "                               capacitated lot-sizing problem with" << endl;
		  out << "                               multiple storage locations (RCLSP-MSL)" << endl;
		  out << endl;
		  out << "Note: In SHP2, the RCLSP-MSL results in the multiple storage location" << endl;
		  out << "problem (MSLP), after the fixing of the decision variables from SHP1 " << endl;
		  out << endl << endl;
          out << "Parameters: Initial Epsilon = " << InitialEpsilon << endl;
		  out << "            Iteration Epsilon = " << Epsilon << endl;
          out << endl << endl << endl;
		  // ****************************************************************************************************



		  

		  //Loop Sequential Heuristic
		  for(;;){


			  //Update one more iteration
			  //to the Sequential Heuristic
			  it_SH += 1;


			  
			  // ************************************ //
			  //    Sequential Heuristic Problem 1    //
			  //             SHP1: RCLSP              //
			  // ************************************ //


			  //Sequential Heuristic Problem 1
			  IloModel SHP1model(env);



			  //Objective Function
			  IloExpr objectiveSHP1(env);


			  //Costs of setup of items
			  for(t=0; t<T; t++)
				 for(i=0; i<I; i++)
					objectiveSHP1 += sc[i]*Y[i][t];


			  //Costs of production
			  for(t=0; t<T; t++)
				 for(tau = t; tau<T; tau++)
					for(i=0; i<I; i++)
					   objectiveSHP1 += vc[i]*FL[i][t][tau];


			  //Cost of inventory of items at storage locations
			  for(t=0; t<T; t++)
				 for(l=0; l<L; l++)
					for(i=0; i<I; i++)
					   objectiveSHP1 += hc[i]*S[i][l][t];


			  //Estimation of the costs of using locations and 
			  //the costs of handling items at storage locations
			  for(t=0; t<T; t++)
				 for(l=0; l<L; l++)
					for(i=0; i<I; i++){
					   objectiveSHP1 += ((g[l]*cs[i])/H[l])*S[i][l][t];

				       objectiveSHP1 += est_ha[i][l]*S[i][l][t];
					} 


			  //Problem objective function environment
			  IloObjective SHP1of = IloMinimize(env, objectiveSHP1);


			  //Add the objective function (SHP1Pof) to the problem
			  SHP1model.add(SHP1of);
			  objectiveSHP1.end();   

			  // *************************************************************




			  //Add the InflowOutflow1 constraints to the problem
			  for(i=0; i<I; i++){
				 InflowOutflow1[i].setNames("InflowOutflow1");
				 SHP1model.add(InflowOutflow1[i]);
			  }





			  //InflowOutflow2SHP1 constraints environment
			  IloArray<IloRangeArray> InflowOutflow2SHP1(env, I);
			  for(i=0; i<I; i++)
				 InflowOutflow2SHP1[i] = IloRangeArray(env, T);

			  //InflowOutflow2SHP1 constraints
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++){
					IloExpr flow(env);

					if (t == 0) {
								   for(tau = t; tau<T; tau++)
									  flow += FL[i][t][tau];

								   flow -= d[i][t];

								   for(l=0; l<L; l++)
									  flow -= S[i][l][t];

								  InflowOutflow2SHP1[i][t] = (flow == 0);
								  flow.end();
					}
					  else {
							  for(tau = t; tau<T; tau++)
								 flow += FL[i][t][tau];

							  flow -= d[i][t];

							  for(l=0; l<L; l++)
								 flow += S[i][l][t-1] - S[i][l][t];

							  InflowOutflow2SHP1[i][t] = (flow == 0);
							  flow.end();
					  }
				 }

			  //Add the InflowOutflow2SHP1 constraints to the problem
			  for(i=0; i<I; i++){
				 InflowOutflow2SHP1[i].setNames("InflowOutflow2SHP1");
				 SHP1model.add(InflowOutflow2SHP1[i]);
			  }





			  //Add the Setup constraints to the problem
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++)
					for(tau = t; tau<T; tau++){
					   Setup[i][t][tau].setName("Setup");
					   SHP1model.add(Setup[i][t][tau]);
					}





			  //Add the Capacity constraints to the problem
			  Capacity.setNames("Capacity");
			  SHP1model.add(Capacity);





			  //CapacityStorageSHP1 constraints environment
    		  IloArray<IloRangeArray> CapacityStorageSHP1(env, L); 
			  for(l=0; l<L; l++){
				 CapacityStorageSHP1[l] = IloRangeArray(env, T);
			  }

			  //CapacityStorageSHP1
			  for(l=0; l<L; l++)
				 for(t=0; t<T; t++){
					IloExpr cap(env);

					for(i=0; i<I; i++)
					   cap += cs[i]*S[i][l][t];
			    
					CapacityStorageSHP1[l][t] = (cap <= H[l]);
					cap.end();
				 }

			  //Add the CapacityStorageSHP1 constraints to the problem
			  for(l=0; l<L; l++){
				 CapacityStorageSHP1[l].setNames("CapacityStorageSHP1");
				 SHP1model.add(CapacityStorageSHP1[l]);
			  }

			  // *************************************************************


			  //Define CPLEX environment to the problem
			  IloCplex SHP1cplex(SHP1model);


			  // *************************************************************



			  //Objective function values, gap and 
		      //computational time found by CPLEX
		      //in the Sequential Heuristic Problem 1 (SHP1)
		      double OF_SHP1 = 0, Gap_SHP1 = 0, Time_SHP1 = 0;


			  //Converte the double variables to binary
			  for(i=0; i<I; i++)
				 SHP1model.add(IloConversion(env, Y[i], ILOBOOL));
	 

		  
			  // ***** Solve the SHP1 ****************************************

			  //Add CPLEX Options 
			  //Print the output and warnings 
			  //of cplex in the output file
			  SHP1cplex.setOut(out);
			  SHP1cplex.setWarning(out);


			  //Limite the number of threads
			  //in the solution of the problem
			  SHP1cplex.setParam(IloCplex::Threads, 1);


			  //Dedecide what CPLEX reports to the screen 
			  //during the solution of the problem
			  SHP1cplex.setParam(IloCplex::MIPDisplay, 3); 

		  		  
			  //Constrol the frequency of displaying node  
			  //logging in the solution of the problem
			  //0 - default: CPLEX's choice; 
			  //n > 0: display new incumbents, and a log line every n nodes
			  SHP1cplex.setParam(IloCplex::MIPInterval, 1000);


			  //Set the maximum time (sec)
			  //Time limite to the SHP1
			  double TimeLimit_SHP1 = (TimeLimit_SH/2);

			  SHP1cplex.setParam(IloCplex::TiLim, TimeLimit_SHP1); 


			  //Set a relative tolerance on the gap between the best 
			  //integer objective and the objective of the best node remaining
			  //CPLEX default: 1e-04 = 0.0001
			  //Pcplex.setParam(IloCplex::EpGap, 0.001); 


	  
			  // ****************************************************************************************************
			  //Print in the output file
			  out << "********************* Solve the SHP1 by an Optimization Package *********************" << endl;
			  out << endl << endl;


			  //SOLVE the problem
			  SHP1cplex.solve(); 

  
			  out << endl << endl << endl;
			  // ****************************************************************************************************



			  // ***** Check the Solution Status *****************************

		  
			  //Check if the solution of the Problem is 
			  //Optimal-1 or Feasible-2
			  if ((SHP1cplex.getStatus() == IloAlgorithm::Optimal) || 
				  (SHP1cplex.getStatus() == IloAlgorithm::Feasible)) {


			  			  //Recover the objective function value
						  OF_SHP1 = SHP1cplex.getValue(SHP1of);


						  //Recover the computational time 
						  //spent to solve the problem
						  Time_SHP1 = SHP1cplex.getTime();


						  //Recover the gap
						  Gap_SHP1 = 100*SHP1cplex.getMIPRelativeGap();



						  //Calculate the costs of SHP1 
						  //without the estimation costs
						  double  CSetupItem_SHP1 = 0, CInventItem_SHP1 = 0;

						  //Cost of setups in the
						  //production of items
						  for(t=0; t<T; t++)
							 for(i=0; i<I; i++)
								if (SHP1cplex.getValue(Y[i][t]) > 0.00001)
								   CSetupItem_SHP1 += sc[i]*(SHP1cplex.getValue(Y[i][t]));

						  //Inventory costs of items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP1cplex.getValue(S[i][l][t]) > 0.00001)
									  CInventItem_SHP1 += hc[i]*(SHP1cplex.getValue(S[i][l][t]));


						  
						  //Calculate the estimation costs in SHP1 
						  double  CSetupLocal_SHP1 = 0, CHandItem_SHP1 = 0;

						  //Estimation costs for using locations 
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP1cplex.getValue(S[i][l][t]) > 0.00001)
								      CSetupLocal_SHP1 += ((g[l]*cs[i])/H[l])*(SHP1cplex.getValue(S[i][l][t]));

						  //Estimation costs for handling items at storage locations
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP1cplex.getValue(S[i][l][t]) > 0.00001)
								      CHandItem_SHP1 += est_ha[i][l]*(SHP1cplex.getValue(S[i][l][t]));



						  // ****************************************************************************************************
						  //Print in the output file
						  out << "******************************** Solution to the SHP1 *******************************" << endl;
						  out << endl << endl;
						  out << "Iteration = " << it_SH << endl;
						  out << "Solution Status SHP1 = " << SHP1cplex.getStatus() << endl;
						  out << "Objective Function Value SHP1 = " << OF_SHP1 << endl; 
						  out << "Gap SHP1 = " << Gap_SHP1 << endl;
						  out << "Time SHP1 = " << Time_SHP1 << endl;
						  out << endl << endl;
						  out << "***** Other Values SHP1 *****" << endl;
						  out << endl;
						  out << "Setup Cost Item SHP1 = " << CSetupItem_SHP1 << endl;
						  out << "Inventory Cost Item SHP1 = " << CInventItem_SHP1 << endl;
						  out << "Estimated Handling Cost Item SHP1 = " << CHandItem_SHP1 << endl;
						  out << "Estimated Setup Cost Location SHP1 = " << CSetupLocal_SHP1 << endl;
						  out << endl;
		                  out << "*************************************************************************************" << endl;
						  // ****************************************************************************************************



			  }//end if Optimal-1/Feasible-2

				else {				

						  // ****************************************************************************************************
						  //Print in the output file
						  out << "************ Upper Bound: The End of the Sequential Heuristic - Version 3 ***********" << endl;
						  out << endl << endl << endl << endl << endl << endl << endl;
						  out << "*************** Final Solution to the Facility Location Reformulation ***************" << endl;
						  out << "***************     of the General Capacitated Lot-Sizing Problem     ***************" << endl;
						  out << "***************      with Multiple Storage Locations (RCLSP-MSL)      ***************" << endl;
						  out << "***************            solved by a Sequential Heuristic           ***************" << endl;
						  out << "***************                       Version 3                       ***************" << endl;
						  out << endl << endl;
						  out << "Solution Status SHP1 = " << SHP1cplex.getStatus() << endl;
						  out << "NO Solution to the RCLSP-MSL" << endl;
						  out << endl << endl;
						  out << "*************************************************************************************" << endl;
						  // ****************************************************************************************************


						  break; //STOP the Sequential Heuristic due to no solution for SHP1

				}//end else


			  // ************************* //
			  //    The End SHP1: RCLSP    //
			  // ************************* //


			  // *************************************************************


			
			  //Create the parameter to recover the
			  //values of the decision variables in 
			  //the optimal solution of the SHP1

			  //Setup of item i in period t
			  IloArray<IloNumArray> Y_sol(env, I);
			  for(i=0; i<I; i++)
				 Y_sol[i] = IloNumArray (env, T);


			  //Inventory of item i at location l in period t
			  IloArray<IloArray<IloNumArray> > S_sol(env, I);
			  for(i=0; i<I; i++){
				 S_sol[i] = IloArray<IloNumArray> (env, L); 
				 for(l=0; l<L; l++){
					S_sol[i][l] = IloNumArray(env, T);
				 }							
			  }


			  //Facility location reformulation 
			  IloArray<IloArray<IloNumArray> > FL_sol(env, I);
			  for(i=0; i<I; i++){
				 FL_sol[i] = IloArray<IloNumArray> (env, T); 
				 for(t=0; t<T; t++){
					FL_sol[i][t] = IloNumArray(env, T);
				 }							
			  }
		  
		  

			  //Initialize the parameters == 0
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++){
					Y_sol[i][t] = 0;

					for(l=0; l<L; l++)
					   S_sol[i][l][t] = 0;

					for(tau = t; tau<T; tau++)
					   FL_sol[i][t][tau] = 0;
				 }
			   


			  //Recover the solution from the SHP1
			  //Setup
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++)
					Y_sol[i][t] = SHP1cplex.getValue(Y[i][t]);


			  //Inventory
			  for(i=0; i<I; i++)
				 for(l=0; l<L; l++)
					for(t=0; t<T; t++)
					   S_sol[i][l][t] = SHP1cplex.getValue(S[i][l][t]);


			  //Reformulation
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++)
					for(tau = t; tau<T; tau++)
					   FL_sol[i][t][tau] = SHP1cplex.getValue(FL[i][t][tau]);

			  // *************************************************************





			  // ************************************ //
			  //    Sequential Heuristic Problem 2    //
			  //            SHP2: RCLSP-MSL           //
			  // ************************************ //


			  //Sequential Heuristic Problem 2
			  IloModel SHP2model(env);



			  //Objective Function
			  IloExpr objectiveSHP2(env);


			  //Costs of setup of items
			  for(t=0; t<T; t++)
				 for(i=0; i<I; i++)
					objectiveSHP2 += sc[i]*Y[i][t];


			  //Costs of production
			  for(t=0; t<T; t++)
				 for(tau = t; tau<T; tau++)
					for(i=0; i<I; i++)
					   objectiveSHP2 += vc[i]*FL[i][t][tau];


			  //Cost of inventory and handling of items at storage locations
			  for(t=0; t<T; t++)
				 for(l=0; l<L; l++)
					for(i=0; i<I; i++)
					   objectiveSHP2 += (hc[i]*S[i][l][t] + ha[i][l]*Dp[i][l][t]);


			  //Costs of using locations
			  for(t=0; t<T; t++)
				 for(l=0; l<L; l++)
					objectiveSHP2 += g[l]*Z[l][t];


			  //Cost of relocation of items between locations
			  for(t=0; t<T; t++)
				 for(l=0; l<L; l++)
					for(k=0; k<L; k++)
					   for(i=0; i<I; i++)
						  objectiveSHP2 += r[i][l][k]*V[i][l][k][t];


			  //Problem objective function environment
			  IloObjective SHP2of = IloMinimize(env, objectiveSHP2);


			  //Add the objective function (SHP2Pof) to the problem
			  SHP2model.add(SHP2of);
			  objectiveSHP2.end();   

			  // *************************************************************



			  //Add the InflowOutflow1 constraints to the problem
			  for(i=0; i<I; i++){
				 InflowOutflow1[i].setNames("InflowOutflow1");
				 SHP2model.add(InflowOutflow1[i]);
			  }





			  //Add the InflowOutflow2 constraints to the problem
			  for(i=0; i<I; i++){
				 InflowOutflow2[i].setNames("InflowOutflow2");
				 SHP2model.add(InflowOutflow2[i]);
			  }





			  //Add the BalanceLocation constraint to the problem
			  for(i=0; i<I; i++)
				 for(l=0; l<L; l++){
					BalanceLocation[i][l].setNames("BalanceLocation");
					SHP2model.add(BalanceLocation[i][l]);
				 }




			 
			  //Add the Setup constraints to the problem
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++)
					for(tau = t; tau<T; tau++){
					   Setup[i][t][tau].setName("Setup");
					   SHP2model.add(Setup[i][t][tau]);
					}





			  //Add the Capacity constraints to the problem
			  Capacity.setNames("Capacity");
			  SHP2model.add(Capacity);



			

			  //Add the InvAlloc constraints to the problem
			  for(i=0; i<I; i++)
				 for(l=0; l<L; l++){
					InvAlloc[i][l].setNames("InvAlloc");
					SHP2model.add(InvAlloc[i][l]);
				 }





			  //Add the CapacityStorage constraints to the problem
			  for(l=0; l<L; l++){
				 CapacityStorage[l].setNames("CapacityStorage");
				 SHP2model.add(CapacityStorage[l]);
			  }





			  //Add the ItemLocatCompat constraints to the problem
			  for(i=0; i<I; i++)
				 for(l=0; l<L; l++){
					ItemLocatCompat[i][l].setNames("ItemLocatCompat");
					SHP2model.add(ItemLocatCompat[i][l]);
				 }





			  //Add the ItemItemCompat constraints to the problem
			  for(i=0; i<I; i++)
				 for(j=i; j<I; j++)
					for(l=0; l<L; l++){
					   ItemItemCompat[i][j][l].setNames("ItemItemCompat");
					   SHP2model.add(ItemItemCompat[i][j][l]);
					}





			  //Add the constraints that fix the decision 
			  //variables in SHP2 according to the values
			  //found in the optimal solution of the SHP1

			  //Setup of item i in period t
			  for(t=0; t<T; t++)
				 for(i=0; i<I; i++)
					SHP2model.add(Y[i][t] == Y_sol[i][t]);



			  //Inventory of item i at location l in period t
			  //RestFixS constraints environment
			  IloArray<IloRangeArray> RestFixS(env, I); 
			  for(i=0; i<I; i++)
				 RestFixS[i] = IloRangeArray(env, T);

			  //RestFixS constraints
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++){
					IloExpr fix(env);
					
					for(l=0; l<L; l++)
					   fix += S[i][l][t];

					for(l=0; l<L; l++)
					   fix -= S_sol[i][l][t];
						   
					RestFixS[i][t] = (fix == 0);
					fix.end();
				 }				

			  //Add the RestFixS constraints to the problem
			  for(i=0; i<I; i++){
				 RestFixS[i].setNames("FixS");
				 SHP2model.add(RestFixS[i]);
			  }



			  //Facility Location Reformulation
			  for(i=0; i<I; i++)
				 for(t=0; t<T; t++)
					for(tau = t; tau<T; tau++)
					   SHP2model.add(FL[i][t][tau] == FL_sol[i][t][tau]);
				

			  // *************************************************************


			  //Define CPLEX environment to the problem
			  IloCplex SHP2cplex(SHP2model);


			  // *************************************************************



			  //Objective function value, gap and 
			  //computational time found by CPLEX
			  //in the Sequential Heuristic Problem 2 (SHP2)
			  double OF_SHP2 = 0, Gap_SHP2 = 0, Time_SHP2 = 0;						

 
			  //Converte the double variables to binary
			  for(l=0; l<L; l++)
				 SHP2model.add(IloConversion(env, Z[l], ILOBOOL));

			  for(i=0; i<I; i++)
				 for(l=0; l<L; l++)
					SHP2model.add(IloConversion(env, W[i][l], ILOBOOL));
		

		  
			  // ***** Solve the SHP2 ****************************************

			  //Add CPLEX Options 
			  //Print the output and warnings 
			  //of cplex in the output file
			  SHP2cplex.setOut(out);
			  SHP2cplex.setWarning(out);


			  //Limite the number of threads
			  //in the solution of the problem
			  SHP2cplex.setParam(IloCplex::Threads, 1);


			  //Dedecide what CPLEX reports to the screen 
			  //during the solution of the problem
			  SHP2cplex.setParam(IloCplex::MIPDisplay, 3); 

		  		  
			  //Constrol the frequency of displaying node  
			  //logging in the solution of the problem
			  //0 - default: CPLEX's choice; 
			  //n > 0: display new incumbents, and a log line every n nodes
			  SHP2cplex.setParam(IloCplex::MIPInterval, 1000);


			  //Set the maximum time (sec)
			  //Time limite to the SHP2
			  double TimeLimit_SHP2 = TimeLimit_SH - Time_SHP1;

			  SHP2cplex.setParam(IloCplex::TiLim, TimeLimit_SHP2); 


			  //Set a relative tolerance on the gap between the best 
			  //integer objective and the objective of the best node remaining
			  //CPLEX default: 1e-04 = 0.0001
			  //Pcplex.setParam(IloCplex::EpGap, 0.001); 


	  
			  // ****************************************************************************************************
			  //Print in the output file
			  out << endl << endl << endl << endl << endl;
			  out << "********************* Solve the SHP2 by an Optimization Package *********************" << endl;
			  out << endl << endl;


			  //SOLVE the problem
			  SHP2cplex.solve(); 

  
			  out << endl << endl << endl;
			  // ****************************************************************************************************



			  // ***** Check the Solution Status *****************************

		  
			  //Check if the solution of the Problem is 
			  //Optimal-1 or Feasible-2
			  if ((SHP2cplex.getStatus() == IloAlgorithm::Optimal) || 
				  (SHP2cplex.getStatus() == IloAlgorithm::Feasible)) {


			  			  //Recover the objective function value
						  OF_SHP2 = SHP2cplex.getValue(SHP2of);


						  //Recover the computational time 
						  //spent to solve the problem
						  Time_SHP2 = SHP2cplex.getTime();


						  //Recover the gap
						  Gap_SHP2 = 100*SHP2cplex.getMIPRelativeGap();


						  
						  //Calculate the costs of SHP1 
						  //without the estimation costs
						  double  CSetupItem_SHP2 = 0, CInventItem_SHP2 = 0, CHandItem_SHP2 = 0, 
							      CSetupLocal_SHP2 = 0, CRelocItem_SHP2 = 0;

						  //Cost of setups of items
						  for(t=0; t<T; t++)
							 for(i=0; i<I; i++)
								if (SHP2cplex.getValue(Y[i][t]) > 0.00001)
								   CSetupItem_SHP2 += sc[i]*(SHP2cplex.getValue(Y[i][t]));

						  //Cost of inventoried items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(S[i][l][t]) > 0.00001)
									  CInventItem_SHP2 += hc[i]*(SHP2cplex.getValue(S[i][l][t]));

						  //Cost of handled items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(Dp[i][l][t]) > 0.00001)
									  CHandItem_SHP2 += ha[i][l]*(SHP2cplex.getValue(Dp[i][l][t]));

						  //Cost of setup used locations
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								if (SHP2cplex.getValue(Z[l][t]) > 0.00001)
								   CSetupLocal_SHP2 += g[l]*(SHP2cplex.getValue(Z[l][t]));

						  //Cost and number of relocated items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(k=0; k<L; k++)
								   for(i=0; i<I; i++)
									  if ((l != k) && (SHP2cplex.getValue(V[i][l][k][t]) > 0.00001))
										 CRelocItem_SHP2 += r[i][l][k]*(SHP2cplex.getValue(V[i][l][k][t]));



						  // ****************************************************************************************************
						  //Print in the output file
						  out << "******************************** Solution to the SHP2 *******************************" << endl;
						  out << endl << endl;
						  out << "Iteration = " << it_SH << endl;
						  out << "Solution Status SHP2 = " << SHP2cplex.getStatus() << endl;
						  out << "Objective Function Value SHP2 = " << OF_SHP2 << endl; 
						  out << "Gap SHP2 = " << Gap_SHP2 << endl;
						  out << "Time SHP2 = " << Time_SHP2 << endl;
						  out << endl << endl;
						  out << "***** Other Values SHP2 *****" << endl;
						  out << endl;
						  out << "Setup Cost Item SHP2 = " << CSetupItem_SHP2 << endl;
						  out << "Inventory Cost Item SHP2 = " << CInventItem_SHP2 << endl;
						  out << "Handling Cost Item SHP2 = " << CHandItem_SHP2 << endl;
						  out << "Setup Cost Location SHP2 = " << CSetupLocal_SHP2 << endl;
						  out << "Relocation Cost Item SHP2 = " << CRelocItem_SHP2 << endl;
						  out << endl;
						  out << "*************************************************************************************" << endl;
						  out << endl << endl << endl << endl << endl;
						  // ****************************************************************************************************


			  }//end if Optimal-1/Feasible-2

				else {				

						  // ****************************************************************************************************
						  //Print in the output file
						  out << "************ Upper Bound: The End of the Sequential Heuristic - Version 3 ***********" << endl;
						  out << endl << endl << endl << endl << endl << endl << endl;
						  out << "*************** Final Solution to the Facility Location Reformulation ***************" << endl;
						  out << "***************     of the General Capacitated Lot-Sizing Problem     ***************" << endl;
						  out << "***************      with Multiple Storage Locations (RCLSP-MSL)      ***************" << endl;
						  out << "***************            solved by a Sequential Heuristic           ***************" << endl;
						  out << "***************                       Version 3                       ***************" << endl;
						  out << endl << endl;
						  out << "Solution Status SHP2 = " << SHP2cplex.getStatus() << endl;
						  out << "NO Solution to the RCLSP-MSL" << endl;
						  out << endl << endl;
						  out << "*************************************************************************************" << endl;
						  // ****************************************************************************************************


						  break; //STOP the Sequential Heuristic due to no solution for SHP2

				}//end else


			  // ***************************** //
			  //    The End SHP2: RCLSP-MSL    //
			  // ***************************** //


			  // *************************************************************





			  //Add the computational time of each  
			  //iteration to the computational 
			  //time of the Sequential Heuristics
			  Time_SH += Time_SHP1 + Time_SHP2;
			  


			  //Save the best solution between
			  //two consecutive iterations 
			  //of the Sequential Heuristic
			  if (SHP2cplex.getValue(SHP2of) <= OF_SH) {


						  //Objective function value
						  OF_SH = SHP2cplex.getValue(SHP2of);


						  //Calculate the gap considering the
						  //linear relaxation as the lower bound
						  Gap_SH = 100*((OF_SH - OF_LR_RCLSPMSL)/OF_SH);



						  // ************************ //
						  //    Other Calculations    // 
						  // ************************ //


						  //Cost and number of setups of items
						  for(t=0; t<T; t++)
							 for(i=0; i<I; i++)
								if (SHP2cplex.getValue(Y[i][t]) > 0.00001) {
								   CSetupItem += sc[i]*(SHP2cplex.getValue(Y[i][t]));

								   NSetupItem += SHP2cplex.getValue(Y[i][t]);
								}


						  //Cost and number of inventoried items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(S[i][l][t]) > 0.00001){
									  CInventItem += hc[i]*(SHP2cplex.getValue(S[i][l][t]));
					  								  
									  NInventItem += SHP2cplex.getValue(S[i][l][t]);
								   }


						  //Cost and number of handled items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(Dp[i][l][t]) > 0.00001) {
									  CHandItem += ha[i][l]*(SHP2cplex.getValue(Dp[i][l][t]));

									  NHandItem += SHP2cplex.getValue(Dp[i][l][t]);
								   }


						  //Cost and number of setup used locations
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								if (SHP2cplex.getValue(Z[l][t]) > 0.00001) {
								   CSetupLocal += g[l]*(SHP2cplex.getValue(Z[l][t]));

								   NLocalUsed += SHP2cplex.getValue(Z[l][t]);
								}


						  //Cost and number of relocated items
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(k=0; k<L; k++)
								   for(i=0; i<I; i++)
									  if ((l != k) && (SHP2cplex.getValue(V[i][l][k][t]) > 0.00001)) {
										 CRelocItem += r[i][l][k]*(SHP2cplex.getValue(V[i][l][k][t]));

										 NRelocItem += SHP2cplex.getValue(V[i][l][k][t]);
									  }

								

						  //Total available space of opened locations
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								TotalOpenSpace += H[l]*(SHP2cplex.getValue(Z[l][t]));
		  

						  //Total used space of the opened locations;
						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   TotalUsedSpace += cs[i]*(SHP2cplex.getValue(S[i][l][t]));
		  
		  
						  //Percentage of used space
						  //over all the opened space
						  PercUsedSpace = 100*(TotalUsedSpace/TotalOpenSpace); 

			  }//end if OF_SH




			  
			  //Objective function value of
			  //the Sequential Heuristic in it2
			  OF_SH_it2 = SHP2cplex.getValue(SHP2of);
			  


			  //STOP the Sequential Heuristic by
			  //the number of iterations or
			  //time limit or
			  //objective function variations
			  if ((it_SH >= 100) || 
				  (TimeLimit_SH < (Time_SHP1 + Time_SHP2)) ||
				  (abs(OF_SH_it1 - OF_SH_it2) < 0.001)) {


						  // ****************************************************************************************************
						  //Print in the output file
						  out << "************ Upper Bound: The End of the Sequential Heuristic - Version 3 ***********" << endl;
						  out << endl << endl << endl << endl << endl << endl << endl;
						  out << "*************** Final Solution to the Facility Location Reformulation ***************" << endl;
						  out << "***************     of the General Capacitated Lot-Sizing Problem     ***************" << endl;
						  out << "***************      with Multiple Storage Locations (RCLSP-MSL)      ***************" << endl;
						  out << "***************            solved by a Sequential Heuristic           ***************" << endl;
						  out << "***************                       Version 3                       ***************" << endl;
						  out << endl << endl;
						  out << "Total Number Iterations = " << it_SH << endl; 
						  out << "Objective Function Value = " << OF_SH << endl; 
						  out << "Gap = " << Gap_SH << endl;
						  out << "Time = " << Time_SH << endl;
						  out << endl << endl;
						  out << "*************************************************************************************" << endl;
						  out << endl << endl << endl;
						  out << "************************************ Other Values ***********************************" << endl;
						  out << endl;
						  out << "Setup Cost Item = " << CSetupItem << endl;
						  out << "Inventory Cost Item = " << CInventItem << endl;
						  out << "Handling Cost Item = " << CHandItem << endl;
						  out << "Setup Cost Location = " << CSetupLocal << endl;
						  out << "Relocation Cost Item = " << CRelocItem << endl;
						  out << endl;
						  out << "Number Setup Item = " << NSetupItem << endl;
						  out << "Number Inventoried Item = " << NInventItem << endl;
						  out << "Number Handled Item = " << NHandItem << endl;
						  out << "Number Used Location = " << NLocalUsed << endl;
						  out << "Number Relocated Item = " << NRelocItem << endl;
						  out << endl;
						  out << "Total Opened Space = " << TotalOpenSpace << endl;
						  out << "Total Used Space = " << TotalUsedSpace << endl;
						  out << "Percentage Used Space = " << PercUsedSpace << endl;
						  out << endl;
						  out << "*************************************************************************************" << endl;
						  out << endl;
						  out << endl;
						  out << endl;
						  out << endl;
						  out << endl;
						  for(t=0; t<T; t++)
							 for(i=0; i<I; i++)
								if (SHP2cplex.getValue(Y[i][t]) > 0.00001)
								   out << "Y_" << i+1 << "_" << t+1 << " = " << SHP2cplex.getValue(Y[i][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(S[i][l][t]) > 0.00001)
									  out << "S_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << SHP2cplex.getValue(S[i][l][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								if (SHP2cplex.getValue(Z[l][t]) > 0.00001)
								   out << "Z_" << l+1 << "_" << t+1 << " = " << SHP2cplex.getValue(Z[l][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(W[i][l][t]) > 0.00001)
									  out << "W_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << SHP2cplex.getValue(W[i][l][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(Dp[i][l][t]) > 0.00001)
									  out << "Dp_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << SHP2cplex.getValue(Dp[i][l][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(Dm[i][l][t]) > 0.00001)
									  out << "Dm_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << SHP2cplex.getValue(Dm[i][l][t]) << endl;

						  out << endl << endl;

						  for(t=0; t<T; t++)
							 for(l=0; l<L; l++)  
								for(k=0; k<L; k++)  
								   for(i=0; i<I; i++)
									  if ((l != k) && (SHP2cplex.getValue(V[i][l][k][t]) > 0.00001))
										 out << "V_" << i+1 << "_" << l+1 << "_" << k+1 << "_" << t+1 << " = " << SHP2cplex.getValue(V[i][l][k][t]) << endl;

						  out << endl << endl;
		  
						  for(t=0; t<T; t++)
							 for(tau = t; tau<T; tau++)
								for(i=0; i<I; i++)
								   if (SHP2cplex.getValue(FL[i][t][tau]) > 0.00001)
									  out << "FL_" << i+1 << "_" << t+1 << "_" << tau+1 << " = " << SHP2cplex.getValue(FL[i][t][tau]) << endl;

						  // ****************************************************************************************************


						  break; //STOP the Sequential Heuristic by Optimality/Feasibility

			  }

			  // *********************************************************
					  




			  //Update the time remaining to solve each 
			  //iteration of the the Sequential Heuristic
			  TimeLimit_SH = TimeLimit_SH - Time_SHP1 - Time_SHP2;

 

			  //Update the objective function value to two
			  //consecutive iterations of the Sequential Heuristic
			  OF_SH_it1 = OF_SH_it2;



			  //Update the estimation of the handling
			  //cost in SHP1 from the solution of the SHP2
			  double NumInvI = 0, CostHandI = 0;			  

			  //Calculate the cost of handling items in inventory
			  for(t=0; t<T; t++)
			     for(l=0; l<L; l++)
				    for(i=0; i<I; i++)
			           if (SHP2cplex.getValue(Dp[i][l][t]) > 0.0001)
				          CostHandI += ha[i][l]*(SHP2cplex.getValue(Dp[i][l][t]));

			  //Calculate the number inventoried items
			  for(t=0; t<T; t++)
			     for(l=0; l<L; l++)
				    for(i=0; i<I; i++)
					   if (SHP1cplex.getValue(S[i][l][t]) > 0.00001)
						  NumInvI += SHP1cplex.getValue(S[i][l][t]);


			  //Update
			  if (NumInvI == 0) {
				  				   for(i=0; i<I; i++)
								      for(l=0; l<L; l++)
				              		     est_ha[i][l] = InitialEpsilon*ha[i][l];
			  }
			    else {				
				        for(l=0; l<L; l++)
				           for(i=0; i<I; i++)
				              est_ha[i][l] = (Epsilon*(CostHandI/NumInvI));				 
				}

			  // *********************************************************



			  //End the parameters of 
			  //the Sequential Heuristic
			  SHP1model.end();
			  SHP1of.end();
			  SHP1cplex.end();
			  SHP2model.end();
			  SHP2of.end();
			  SHP2cplex.end();



		  }//end Loop for ;;
		  //Sequential Heuristic


		  // ***************************************** //
		  //    The End of the Sequential Heuristic    //
		  // ***************************************** //


					  


// ********************************************************************************************** //
// *************************************** END MAIN PROGRAM ************************************* //
// ********************************************************************************************** //



   } //end try


      catch (IloException& ex) {
        cerr << "Error Cplex: " << ex << endl;
      } 
	  catch (...) {
        cerr << "Error Cpp" << endl;
      }
     

	  //Finish the enviroment
	  env.end();

	
	return 0;

}

