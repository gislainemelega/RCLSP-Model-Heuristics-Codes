// ************************************************************************************* //
//    Program to solve the Facility Location Reformulation of the General Capacitated    //
//             Lot-Sizing Problem with Multiple Storage Locations (RCLSP-MSL)            //
//								by an Optimization Package								 //
//	 																					 //
//	  Used in Gislaine Mara Melega Pos-doctoral											 //
//    Copyright 2020 - date:  04/2022													 //
// ************************************************************************************* //



//Libraries
#include <stdafx.h>
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


int main(int, char **)
{

	//Input Data File name
    ifstream in("dataRCLSPMSL.dat");


	//Output Data File name
	ofstream out("resultsRCLSPMSL.dat");


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



// *********************************** // 
//      Creation of the Parameters     //
//    and Data Read from Input File    // 
// *********************************** //


		  //Indexes
		  //Time period
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
		  //Production cost per unit of item
          IloNumArray vc(env, I);

          //Setup cost of item 
		  IloNumArray sc(env, I);										

		  //Inventory holding cost per unit of item
		  IloNumArray hc(env, I);

		  //Production consumption of capacity per unit of item 
		  IloNumArray vt(env, I);

		  //Setup time consumption of capacity of item 
		  //IloNumArray st(env, I);

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

		  //Production of item i in period t
		  IloArray<IloNumVarArray> X(env, I);
		  for(i=0; i<I; i++){
		     X[i] = IloNumVarArray(env, T);
			 for(t=0; t<T; t++){
                sprintf_s(proditem, "X_%d_%d", i, t);
				X[i][t] = IloNumVar(env, 0, IloInfinity, proditem);
			 }
		  }



		  //Setup of item i in period t
		  IloArray<IloNumVarArray> Y(env, I);
		  for(i=0; i<I; i++){
		     Y[i] = IloNumVarArray(env, T);
			 for(t=0; t<T; t++){
                sprintf_s(setupitem, "Y_%d_%d", i, t);
				Y[i][t] = IloNumVar(env, 0, 1, setupitem);
			 }
		  }



		  //Inventory of item i at location l in period t
		  IloArray<IloArray<IloNumVarArray> > S(env, I);
		  for(i=0; i<I; i++){
		     S[i] = IloArray<IloNumVarArray> (env, L); 
		     for(l=0; l<L; l++){ 
		        S[i][l] = IloNumVarArray(env, T);								
				for(t=0; t<T; t++){ 
				   sprintf_s(stockitemloc, "S_%d_%d_%d", i, l, t);
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
                sprintf_s(setuplocation, "Z_%d_%d", l, t);
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
				   sprintf_s(inflowitem, "Dp_%d_%d_%d", i, l, t);
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
				   sprintf_s(outflowitem, "Dm_%d_%d_%d", i, l, t);
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
				   sprintf_s(assignitemlocation, "W_%d_%d_%d", i, l, t);
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
					   sprintf_s(relocationitem, "V_%d_%d_%d_%d", i, l, k, t);
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
				   sprintf_s(facilityref, "FL_%d_%d_%d", i, t, tau);
                   FL[i][t][tau] = IloNumVar(env, 0, IloInfinity, facilityref);
				}
			 }
		  }

		  // *************************************************************



		  //Objective Function
		  IloExpr objective(env);


		  //Costs of production
		  for(t=0; t<T; t++)
			 for(tau = t; tau<T; tau++)
			    for(i=0; i<I; i++)
		           objective += vc[i]*FL[i][t][tau];


		  //Costs of setup of items
		  for(t=0; t<T; t++)
			 for(i=0; i<I; i++)
		        objective += sc[i]*Y[i][t];


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



		  //InflowOutFlow1 constraints environment
		  IloArray<IloRangeArray> InflowOutFlow1(env, I);
		  for(i=0; i<I; i++)
		     InflowOutFlow1[i] = IloRangeArray(env, T);

		  //InflowOutFlow1 constraints
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++){
				IloExpr flow(env);
				
				for(tau=0; tau <= t; tau++)
				   flow += FL[i][tau][t];

			    InflowOutFlow1[i][t] = (flow == d[i][t]);
				flow.end();
			 }

		  //Add the InflowOutFlow1 constraints to the problem
		  for(i=0; i<I; i++){
             InflowOutFlow1[i].setNames("InflowOutFlow1");
	         Pmodel.add(InflowOutFlow1[i]);
		  }





		  //InflowOutFlow2 constraints environment
		  IloArray<IloRangeArray> InflowOutFlow2(env, I);
		  for(i=0; i<I; i++)
		     InflowOutFlow2[i] = IloRangeArray(env, T);

		  //InflowOutFlow2 constraints
		  for(i=0; i<I; i++)
		     for(t=0; t<T; t++){
				IloExpr flow(env);
				
				for(tau = t; tau<T; tau++)
				   flow += FL[i][t][tau];

				flow -= d[i][t];

				for(l=0; l<L; l++)
				   flow -= Dp[i][l][t] - Dm[i][l][t];

			    InflowOutFlow2[i][t] = (flow == 0);
				flow.end();
			 }

		  //Add the InflowOutFlow2 constraints to the problem
		  for(i=0; i<I; i++){
             InflowOutFlow2[i].setNames("InflowOutFlow2");
	         Pmodel.add(InflowOutFlow2[i]);
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
				   double allocitem = 0;

				   //Sum of the demand
				   for(tau=t; tau<T; tau++)
				      sumd += d[i][tau];

				   //Allocation of items
				   allocitem = (H[l]/cs[i]);

				   //Choose the smallest value between the sum
				   //of the demand and allocation of items
				   if (sumd <= allocitem) BigM = sumd;
				     else BigM = allocitem;


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
		  //the solution of the linear relation
		  Pcplex.setParam(IloCplex::Threads, 1);



		  //Extract the model
		  //Pcplex.extract(Pmodel);
		  //Export the model
		  //Pcplex.exportModel("model.lp" );

	 
		  
          // ****************************************************************************************************
		  out << "***** The Facility Location Reformulation of the General Capacitated Lot-Sizing *****" << endl;
		  out << "*****            Problem with Multiple Storage Locations (RCLSP-MSL)            *****" << endl;       
		  out << "*****                     solved by an Optimization Package                     *****" << endl;
		  out << endl << endl;
		  out << "* Lower Bound: Linear Relaxation" << endl;
		  out << "* Upper Bound: Optimization Package" << endl;
		  out << endl << endl << endl;
		  out << "*********** Lower Bound: Solve the LR_RCLSP-MSL by an Optimization Package **********" << endl;
		  out << endl << endl;


		  //SOLVE the LR_RCLSP-MSL
          Pcplex.solve();


		  out << endl << endl << endl << endl << endl;
		  // ****************************************************************************************************	  





		  // ***** Check the Solution Status *****************************

		  
		  //Check if the solution of the Linear Relaxation is
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
					out << endl << endl << endl << endl << endl;
					// ****************************************************************************************************


		  }//end if Optimal-1/Feasible-2
		    
		    else {

				    // ****************************************************************************************************
					//Print in the output file
					out << "*********** Solution to the Linear Relaxation of the General Capacitated ************" << endl;
					out << "***********      Lot-Sizing Problem with Multiple Storage Locations      ************" << endl;
					out << endl << endl;
					out << "Solution Status LR = " << Pcplex.getStatus() << endl;
					out << "NO Solution to the Linear Relaxation of the RCLSP-MSL" << endl;
					out << endl << endl;
					out << "*************************************************************************************" << endl;
					out << endl << endl << endl << endl << endl << endl << endl << endl << endl << endl;
					// ****************************************************************************************************

			}//end else

// ***********************************************************************





// ******************************************************************* //
//    Solve the Reformulated General Capacitated Lot-Sizing Problem    // 
//			   with Multiple Storage Locations (RCLSP-MSL)			   //
//                     by an Optimization Package                      //
// ******************************************************************* //


		  //Objective function values, gap and 
		  //computational time found by CPLEX
		  double OF_BestLB_RCLSPMSL, OF_RCLSPMSL, Gap_RCLSPMSL, Time_RCLSPMSL;									


		  //Converte the double variables to binary
		  for(i=0; i<I; i++)
		     Pmodel.add(IloConversion(env, Y[i], ILOBOOL));
	 
		  for(l=0; l<L; l++)
			 Pmodel.add(IloConversion(env, Z[l], ILOBOOL));

		  for(i=0; i<I; i++)
		     for(l=0; l<L; l++)
		        Pmodel.add(IloConversion(env, W[i][l], ILOBOOL));
		

		  
		  // ***** Solve the RCLSP-MSL ***********************************

		  //Add CPLEX Options 
          //Print the output and warnings 
		  //of cplex in the output file
          Pcplex.setOut(out);
          Pcplex.setWarning(out);


		  //Limite the number of threads
		  //in the solution of the problem
		  Pcplex.setParam(IloCplex::Threads, 1);


		  //Dedecide what CPLEX reports to the screen 
		  //during the solution of the problem
		  Pcplex.setParam(IloCplex::MIPDisplay, 3); 

		  		  
		  //Constrol the frequency of displaying node  
		  //logging in the solution of the problem
		  //0 - default: CPLEX's choice; 
		  //n > 0: display new incumbents, and a log line every n nodes
		  Pcplex.setParam(IloCplex::MIPInterval, 1000);


		  //Set the maximum time (seconds) 
		  //to solve the problem
		  //CPLEX default: 1e+75
		  Pcplex.setParam(IloCplex::TiLim, 1800); 


		  //Set a relative tolerance on the gap between the best 
		  //integer objective and the objective of the best node remaining
		  //CPLEX default: 1e-04 = 0.0001
		  //Pcplex.setParam(IloCplex::EpGap, 0.0001); 


	  
		  // ****************************************************************************************************
		  //Print in the output file
		  out << "************ Upper Bound: Solve the RCLSP-MSL by an Optimization Package ************" << endl;
		  out << endl << endl;


		  //SOLVE the RCLSP-MSL
	      Pcplex.solve(); 

  
		  out << endl << endl;
		  out << "*************************************************************************************" << endl;
		  // ****************************************************************************************************





		  // ***** Check the Solution Status *****************************

		  
		  //Check if the solution of the Problem is 
		  //Optimal-1 or Feasible-2
		  if ((Pcplex.getStatus() == IloAlgorithm::Optimal) || 
			  (Pcplex.getStatus() == IloAlgorithm::Feasible)) {


					  //Recover the best lower bound 
				      //found by CPLEX to the SLP
					  OF_BestLB_RCLSPMSL = Pcplex.getBestObjValue();


			  		  //Recover the objective function value
					  OF_RCLSPMSL = Pcplex.getValue(Pof);


					  //Recover the gap
					  Gap_RCLSPMSL = 100*Pcplex.getMIPRelativeGap();


					  //Recover the computational time 
					  //spent to solve the problem
					  Time_RCLSPMSL = Pcplex.getTime();
					 


					  // ************************ //
					  //    Other Calculations    // 
					  // ************************ //


					  double CSetupItem = 0,
						     CInventItem = 0,
							 CHandItem = 0, 
							 CSetupLocal = 0, 
							 CRelocItem = 0,
							 
							 NSetupItem = 0,
							 NInventItem = 0,
							 NHandItem = 0,
							 NLocalUsed = 0,
							 NRelocItem = 0,

							 TotalOpenSpace = 0,
					         TotalUsedSpace = 0,
					         PercUsedSpace = 0;



					  //Cost and number of setups 
					  //in the production of items
					  for(t=0; t<T; t++)
					     for(i=0; i<I; i++)
							if (Pcplex.getValue(Y[i][t]) > 0.00001) {
						       CSetupItem += sc[i]*(Pcplex.getValue(Y[i][t]));

							   NSetupItem += Pcplex.getValue(Y[i][t]);
							}


					  //Inventory costs of items
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
							for(i=0; i<I; i++)
							   if (Pcplex.getValue(S[i][l][t]) > 0.00001){
							      CInventItem += hc[i]*(Pcplex.getValue(S[i][l][t]));

								  NInventItem += Pcplex.getValue(S[i][l][t]);
							   }


					  //Number and cost of handled items
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
							for(i=0; i<I; i++)
							   if (Pcplex.getValue(Dp[i][l][t]) > 0.00001) {
					              CHandItem += ha[i][l]*(Pcplex.getValue(Dp[i][l][t]));

								  NHandItem += Pcplex.getValue(Dp[i][l][t]);
							   }



					  //Number and setup cost of used locations
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
							if (Pcplex.getValue(Z[l][t]) > 0.00001) {
							   CSetupLocal += g[l]*(Pcplex.getValue(Z[l][t]));

							   NLocalUsed += Pcplex.getValue(Z[l][t]);
							}


					  //Number and cost of relocated items
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
							for(k=0; k<L; k++)
							   for(i=0; i<I; i++)
								  if ((l != k) && (Pcplex.getValue(V[i][l][k][t]) > 0.00001)) {
								     CRelocItem += r[i][l][k]*(Pcplex.getValue(V[i][l][k][t]));

								     NRelocItem += Pcplex.getValue(V[i][l][k][t]);
								  }

								  
								  
					  //Total available space of used locations
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        TotalOpenSpace += H[l]*(Pcplex.getValue(Z[l][t]));


					  //Total used space in the locations;
					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        for(i=0; i<I; i++)
						       TotalUsedSpace += cs[i]*(Pcplex.getValue(S[i][l][t]));


					  //Percentage of used space
					  //over all the opened space
					  PercUsedSpace = 100*(TotalUsedSpace/TotalOpenSpace);

	  



					  // ****************************************************************************************************
					  //Print in the output file
					  out << endl << endl << endl << endl << endl;
					  out << "*************** Final Solution to the Facility Location Reformulation ***************" << endl;
					  out << "***************     of the General Capacitated Lot-Sizing Problem     ***************" << endl;
					  out << "***************      with Multiple Storage Locations (RCLSP-MSL)      ***************" << endl;
					  out << endl << endl;
					  out << "Solution Status = " << Pcplex.getStatus() << endl;
					  out << "Number of Nodes = " << Pcplex.getNnodes() << endl;
					  out << "Best Lower Bound Value = " << OF_BestLB_RCLSPMSL << endl;
					  out << "Objective Function Value = " << OF_RCLSPMSL << endl; 
					  out << "Gap = " << Gap_RCLSPMSL << endl;
					  out << "Time = " << Time_RCLSPMSL << endl;
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
						    if (Pcplex.getValue(Y[i][t]) > 0.00001)
							   out << "Y_" << i+1 << "_" << t+1 << " = " << Pcplex.getValue(Y[i][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        for(i=0; i<I; i++)
						       if (Pcplex.getValue(S[i][l][t]) > 0.00001)
							      out << "S_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << Pcplex.getValue(S[i][l][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
					     for(l=0; l<L; l++)
						    if (Pcplex.getValue(Z[l][t]) > 0.00001)
							   out << "Z_" << l+1 << "_" << t+1 << " = " << Pcplex.getValue(Z[l][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        for(i=0; i<I; i++)
						       if (Pcplex.getValue(W[i][l][t]) > 0.00001)
							      out << "W_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << Pcplex.getValue(W[i][l][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        for(i=0; i<I; i++)
						       if (Pcplex.getValue(Dp[i][l][t]) > 0.00001)
							      out << "Dp_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << Pcplex.getValue(Dp[i][l][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)
					        for(i=0; i<I; i++)
						       if (Pcplex.getValue(Dm[i][l][t]) > 0.00001)
							      out << "Dm_" << i+1 << "_" << l+1 << "_" << t+1 << " = " << Pcplex.getValue(Dm[i][l][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
						 for(l=0; l<L; l++)  
						    for(k=0; k<L; k++)  
							   for(i=0; i<I; i++)
						          if ((l != k) && (Pcplex.getValue(V[i][l][k][t]) > 0.00001))
							         out << "V_" << i+1 << "_" << l+1 << "_" << k+1 << "_" << t+1 << " = " << Pcplex.getValue(V[i][l][k][t]) << endl;

					  out << endl << endl;

					  for(t=0; t<T; t++)
					     for(tau = t; tau<T; tau++)
							for(i=0; i<I; i++)
							   if (Pcplex.getValue(FL[i][t][tau]) > 0.00001)
							       out << "FL_" << i+1 << "_" << t+1 << "_" << tau+1 << " = " << Pcplex.getValue(FL[i][t][tau]) << endl;

				      // ******************************************************************************************************************






		  }//end if Optimal-1/Feasible-2
		    
		    //else1
		    else {

					  // ****************************************************************************************************
					  //Print in the output file
					  out << endl << endl << endl << endl << endl;
					  out << "*************** Final Solution to the Facility Location Reformulation ***************" << endl;
					  out << "***************     of the General Capacitated Lot-Sizing Problem     ***************" << endl;
					  out << "***************      with Multiple Storage Locations (RCLSP-MSL)      ***************" << endl;
					  out << endl << endl;
					  out << "Solution Status = " << Pcplex.getStatus() << endl;
					  out << "NO Solution to the RCLSP-MSL" << endl;
					  out << endl << endl;
					  out << "*************************************************************************************" << endl;
					  // ****************************************************************************************************

			}//end else



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
