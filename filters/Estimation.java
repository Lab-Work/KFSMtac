package filters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import triggerer.*;
import model.*;

import org.jblas.DoubleMatrix;

import doubleMatrix.InverseMatrix;
import trueSolution.TrueSolution;
import trueSolution.TrueinGround;
import targetGroundTruth.TargetGroundTruth;

public class Estimation {
	
	Filter filter;

	public Estimation(String _nameFilter, Triggerer _triggerer, ModelKF _modelKF) {
		Filter f = Filter.createFilter(_nameFilter, _triggerer, _modelKF);	
		filter=f;
	}
		
	public void updateTrueSolution() {
		filter.modelKF.updateTrueSolution();
	}	

	public static void exportResult(TargetGroundTruth _targetGroundTruth, int limite, String folder) {
		try {
			
			int numAvg=100;
		    int numKFSM=numAvg;
		    int numKFDT=numAvg;
		    int numKFST=numAvg;
			int numFilters = numKFSM+numKFDT+numKFST;
			int cells=_targetGroundTruth.cells;	
			
		    TrueSolution[] trueSolution=new TrueSolution[numFilters];
		    for(int i=0;i<numFilters;i++){
		        trueSolution[i] = new TrueinGround(_targetGroundTruth,cells,i, numFilters);
		    }
		    
			ModelKF [] modelKFs=new ModelKF[numFilters];
			for (int i=0;i<numFilters;i++){					
				modelKFs[i]=trueSolution[i].setModels();
			}
					
			BufferedWriter writerTrue;
			BufferedWriter writerRate;
			BufferedWriter writerTrVar;
			BufferedWriter[] writerFilter = new BufferedWriter[3];
			BufferedWriter[] writerError= new BufferedWriter[3];
			
			Triggerer []DTtriggerSM=new Triggerer [numKFSM];
			Triggerer []DTtriggerDT=new Triggerer [numKFDT];
			Triggerer []STtriggerST=new Triggerer [numKFST];
			for(int i=0; i<numKFSM;i++){
				DTtriggerSM[i]=Triggerer.createTriggerer("DT");
				DTtriggerSM[i].setNewParameters(DoubleMatrix.eye(1).mul(10),0.00);
			}
			for(int i=0; i<numKFDT;i++){
				DTtriggerDT[i]=Triggerer.createTriggerer("DT");
				DTtriggerDT[i].setNewParameters(DoubleMatrix.eye(1).mul(200000), 0.00);
			}
			for(int i=0; i<numKFST;i++){
				STtriggerST[i]=Triggerer.createTriggerer("ST");
				STtriggerST[i].setNewParameters(DoubleMatrix.eye(1).mul(25000000), 0.31);
			}
						
			System.out.print("Beginning of simulation");
				
				new File("results/"+folder).mkdir();
				
				Estimation[] estimations = new Estimation[numFilters];
				double [] rate =new double [3];
				double [] erroravg =new double [3];
				double [] error1 =new double [numFilters];
				double [][] errorSDV =new double [3][limite];
				double [] errorSDV_avg =new double [3];
				double [] errorSDVk =new double [3];
				double [] errorplusSDV =new double [3];
				int kinitial=30;
				
				writerTrue = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"trueState.csv")));
				writerRate = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Rate.csv")));
				writerTrVar = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Trace.csv")));
				
				for (int i = 0; i<numFilters; i++) {
					if (i<numKFSM){
						estimations[i] = new Estimation("KFSM", DTtriggerSM[i], modelKFs[i]);
					}
					else if (i>=numKFSM && i<numKFSM+numKFDT){
						estimations[i] = new Estimation("KFDT", DTtriggerDT[i-numKFSM], modelKFs[i]);
					}
					else{
						estimations[i] = new Estimation("KFST", STtriggerST[i-numKFSM-numKFDT], modelKFs[i]);
					}

				}

				for (int i = 0; i<3; i++) {
					
					String S = "results/"+folder+"/"+estimations[i*numAvg].filter.getClass().getSimpleName();
					writerFilter[i] = new BufferedWriter(new FileWriter(new File(S+".csv")));
					writerError[i]=new BufferedWriter(new FileWriter(new File(S+"error.csv")));		
				}     
				
				for (int k=0; k<limite; k++) {
					
					for (int i=0;i<_targetGroundTruth.cells;i++){
						writerTrue.write(_targetGroundTruth.trueStatesGround.get(i,0)+",");
					}
					writerTrue.write("\n");		

					DoubleMatrix[] mean = new DoubleMatrix[3];	

					for (int i = 0; i<3; i++) {
						mean[i]= DoubleMatrix.zeros(_targetGroundTruth.trueStatesGround.getRows(), 1);
					}

					for (int i=0; i<3; i++) {
						errorSDVk[i]=0;
						if (i==0){
							double trace=0;
							for (int j=0;j<numKFSM; j++){
								mean[i]=mean[i].add((estimations[j].filter.mean).mul(((double)(1)/(double)(numKFSM))));
								for (int l=0; l<mean[i].length ; l++) {
            						trace=trace+estimations[j].filter.var.get(l,l);
            					}  	
							}
							trace=trace*((double)(1)/(double)(numKFSM));
							writerTrVar.write(trace+",");
							for (int j=0;j<numKFSM; j++){
								errorSDVk[i]=errorSDVk[i]+((((estimations[j].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j].filter.mean.sub(mean[i]))).get(0,0));
							}
							errorSDVk[i]=Math.sqrt(errorSDVk[i])/((double)(numKFSM));
						}
						else if (i==1){
							double trace1=0;
							for (int j=0;j<numKFDT; j++){
								mean[i]=mean[i].add((estimations[j+numKFSM].filter.mean).mul(((double)(1)/(double)(numKFDT))));
								for (int l=0; l<mean[i].length ; l++) {
            						trace1=trace1+estimations[j+numKFSM].filter.var.get(l,l);
            					}  	
							}
							trace1=trace1*((double)(1)/(double)(numKFDT));
							writerTrVar.write(trace1+",");
							for (int j=0;j<numKFDT; j++){
								errorSDVk[i]=errorSDVk[i]+((((estimations[j+numKFSM].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j+numKFSM].filter.mean.sub(mean[i]))).get(0,0));
							}
							errorSDVk[i]=Math.sqrt(errorSDVk[i])/((double)(numKFDT));
						}
						else{
							double trace2=0;
							for (int j=0;j<numKFST; j++){
								mean[i]=mean[i].add((estimations[j+numKFSM+numKFDT].filter.mean).mul(((double)(1)/(double)(numKFST))));
								for (int l=0; l<mean[i].length ; l++) {
            						trace2=trace2+estimations[j+numKFSM+numKFDT].filter.var.get(l,l);
            					}  	
							}
							trace2=trace2*((double)(1)/(double)(numKFST));
							writerTrVar.write(trace2+",");
							for (int j=0;j<numKFST; j++){
								errorSDVk[i]=errorSDVk[i]+((((estimations[j+numKFSM+numKFDT].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j+numKFSM+numKFDT].filter.mean.sub(mean[i])).get(0,0)));
							}
							errorSDVk[i]=Math.sqrt(errorSDVk[i])/((double)(numKFST));
							
						}
																					
						for (int j=0; j<mean[i].length ; j++) {
							writerFilter[i].write(mean[i].get(j)+",");							
						}
						
    					for (int j=0; j<mean[i].length ; j++) {
							writerError[i].write(mean[i].get(j)-_targetGroundTruth.trueStatesGround.get(j)+",");						
						}
    					writerError[i].write(Math.sqrt((((mean[i].sub(_targetGroundTruth.trueStatesGround)).transpose()).mmul(mean[i].sub(_targetGroundTruth.trueStatesGround))).get(0,0))+",");
    					writerError[i].write(errorSDVk[i]+",");

						writerError[i].write("\n");
						writerFilter[i].write("\n");
						if(k>=kinitial){
							erroravg[i]=erroravg[i]+(Math.sqrt((((mean[i].sub(_targetGroundTruth.trueStatesGround)).transpose()).mmul(mean[i].sub(_targetGroundTruth.trueStatesGround))).get(0,0)))/((double)(limite-kinitial));
						}
						
						
					}
					
					writerTrVar.write("\n");
					
					
					if(k>=kinitial){

						for (int i=0; i<numFilters; i++) {
		
							error1[i]=error1[i]+(Math.sqrt((((estimations[i].filter.mean.sub(_targetGroundTruth.trueStatesGround)).transpose()).mmul(estimations[i].filter.mean.sub(_targetGroundTruth.trueStatesGround))).get(0,0)))/((double)(limite-kinitial));
							 
							
						}
						for (int i=0;i<3;i++){
							if (i==0){
								for (int j=0;j<numKFSM; j++){
									errorSDV[i][k]=errorSDV[i][k]+((((estimations[j].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j].filter.mean.sub(mean[i]))).get(0,0));
								}
								errorSDV[i][k]=Math.sqrt(errorSDV[i][k]/((double)(numKFSM)));
							}
							else if (i==1){
								for (int j=0;j<numKFDT; j++){
									errorSDV[i][k]=errorSDV[i][k]+((((estimations[j+numKFSM].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j+numKFSM].filter.mean.sub(mean[i]))).get(0,0));
								}
								errorSDV[i][k]=Math.sqrt(errorSDV[i][k]/((double)(numKFDT)));
							}
							else{
								for (int j=0;j<numKFST; j++){
									errorSDV[i][k]=errorSDV[i][k]+((((estimations[j+numKFSM+numKFDT].filter.mean.sub(mean[i])).transpose()).mmul(estimations[j+numKFSM+numKFDT].filter.mean.sub(mean[i])).get(0,0)));
								}
								errorSDV[i][k]=Math.sqrt(errorSDV[i][k]/((double)(numKFST)));
							}
						}
					}
					
					
					double [] ratesep =new double [numFilters];
					for (int i=0; i<numFilters; i++) {
						ratesep[i]=ratesep[i]+((double)(estimations[i].filter.triggerer.gamma)/(double)(limite));
						if(numAvg==1){
	                    	writerRate.write(estimations[i].filter.triggerer.gamma+",");
	                    }
					
					}
					writerRate.write("\n");

					
					if (k==2){
						estimations[1].filter.var.print();
						estimations[1].filter.f_var.print();
						estimations[1].filter.f_mean.print();
						DoubleMatrix compute=DoubleMatrix.zeros(2,1);
						DoubleMatrix computeK=DoubleMatrix.zeros(2,1);
						DoubleMatrix compute1=DoubleMatrix.zeros(1,1);
						
						compute1=InverseMatrix.invPoSym( (estimations[1].filter.measure.mmul(estimations[1].filter.f_var.mmul(estimations[1].filter.measure.transpose()))).add(estimations[1].filter.measureVar) );
						computeK=estimations[1].filter.f_var.mmul(((estimations[1].filter.measure).transpose()).mmul(compute1));
						computeK.print();
						compute=estimations[1].filter.f_mean.add(computeK.mmul(estimations[1].filter.measurements.sub(estimations[1].filter.measure.mmul(estimations[1].filter.f_mean))));
						mean[1].print();
						compute.print();
						mean[1].sub(_targetGroundTruth.trueStatesGround).print();
						estimations[1].filter.measure.print();
					}
					
					

					
					
					for (int i=0; i<3; i++) {
						if (i==0){
							for (int j=0;j<numKFSM; j++){
								rate[i]=rate[i]+((((double)(ratesep[j])/(double)(numKFSM))));
							}
						}
						else if (i==1){
							for (int j=0;j<numKFDT; j++){
								rate[i]=rate[i]+((((double)(ratesep[j+numKFSM])/(double)(numKFDT))));
							}							
						}
						else{
							for (int j=0;j<numKFST; j++){
								rate[i]=rate[i]+((((double)(ratesep[j+numKFSM+numKFDT])/(double)(numKFST))));
							}	
						}
					}
													
					for (int i=0; i<numFilters; i++) {
						estimations[i].filter.getNewParametersFromModel();
					}			
				
					_targetGroundTruth.update();			

					for(int i=0; i<numFilters; i++){
						trueSolution[i].update();
					}
					
					for(int i=0; i<numFilters; i++){
						trueSolution[i].updatemeasurement();
					}

					for (int i=0; i<numFilters; i++) {
						estimations[i].filter.NextStep();
					}					
		
					if (k==limite-1){
						double error_0=0;
						double error_1=0;
						double error_2=0;
						


						
						for (int i=0; i<3; i++) {
							if (i==0){						
								for (int j=0;j<numKFSM; j++){
									error_0=error_0+error1[j]*((double)(1)/(double)(numKFSM));
								}
							}
							else if (i==1){
								
								for (int j=0;j<numKFDT; j++){
									error_1=error_1+error1[j+numKFSM]*((double)(1)/(double)(numKFDT));
								}
							}
							else{
								
								for (int j=0;j<numKFST; j++){
									error_2=error_2+error1[j+numKFSM+numKFDT]*((double)(1)/(double)(numKFST));
								}
							}
							
							for (int k1=kinitial;k1<=limite-1;k1++){
								errorSDV_avg[i]=errorSDV_avg[i]+(errorSDV[i][k1]/(double)(limite-kinitial));
							}
							
							errorplusSDV[i]=erroravg[i]+errorSDV_avg[i];
						}
						
						
						System.out.print("error_NOcancel_KFSM="+error_0+"\n");
						System.out.print("error_NOcancel_KFDT="+error_1+"\n");
						System.out.print("error_NOcancel_KFST="+error_2+"\n");
						System.out.print("error_avg_KFSM="+erroravg[0]+"\n");
						System.out.print("error_avg_KFDT="+erroravg[1]+"\n");
						System.out.print("error_avg_KFST="+erroravg[2]+"\n");
						System.out.print("error_plus_SDV_KFSM="+errorplusSDV[0]+"\n");
						System.out.print("error_plus_SDV_KFDT="+errorplusSDV[1]+"\n");
						System.out.print("error_plus_SDV_KFST="+errorplusSDV[2]+"\n");
						writerRate.write("rate_KFSM="+rate[0]+"\n");
						writerRate.write("rate_KFDT="+rate[1]+"\n");
						writerRate.write("rate_KFST="+rate[2]+"\n");
					}
				}
				
				for (int i = 0; i<3; i++) {
					writerFilter[i].flush(); writerFilter[i].close();
					writerError[i].flush(); writerError[i].close();
				}

				writerTrue.flush(); writerTrue.close();
				writerRate.flush(); writerRate.close();
				writerTrVar.flush(); writerTrVar.close();
				System.out.print("rate_KFSM="+rate[0]+"\n");
				System.out.print("rate_KFDT="+rate[1]+"\n");
				System.out.print("rate_KFST="+rate[2]+"\n");
			System.out.println(" End");			
		}
		catch (Exception e) {e.printStackTrace();}
		
	}
	
}
					