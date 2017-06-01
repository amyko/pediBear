package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import mcmcMoves.Move;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import utility.DataParser;

public class SimulatedAnnealing {
	
	final Pedigree ped;
	final int runLength;
	final Move[] moves; 
	final PrintWriter cranefootFamWriter;
	//final PrintWriter pairWriter;
	PrintWriter convWriter;
	final Random rGen;
	private int missingParentCounter = 0;
	
	//heating parameters
	private double[] heat;
	private int[] coolingSchedule;
	
	
	
	public double bestLkhd = Double.NEGATIVE_INFINITY;
	public double stopThresh;;
	

	public SimulatedAnnealing(Pedigree ped, double[] heat, int[] coolingSchedule, Move[] moves, int runLength, Random rGen, String outPath, double stopThresh) throws IOException{

		this.ped = ped;
		this.runLength = runLength;
		this.moves = moves;		
		this.cranefootFamWriter = DataParser.openWriter(outPath+".fam");
		//this.pairWriter = DataParser.openWriter(outPath+".pair");
		this.convWriter = DataParser.openWriter(outPath+".lkhd");
		
		this.rGen = rGen;
		this.heat = heat;
		this.coolingSchedule = coolingSchedule;
		this.stopThresh = stopThresh;
		
	}
	
	
	
	public void run(){

		runBurnIn();
		
		runSample();
	
	}

	
	
	
	private void runBurnIn(){
		
		//System.out.println("Running simulated annealing");
		
		convWriter.write(">\n");

		
		int iter = 0;
		double initLkhd = 0;
		double counter = 0;
		
		//for each temperature
		for(int t=0; t<heat.length-1; t++){
			
			if(iter%100000==0){
				initLkhd = ped.getLogLikelihood();
				counter = 0;
			}

			//run burn-in
			//while(moves[0].nTried < 100 && moves[0].nAccept <50){
			for(int i = 0; i < coolingSchedule[t]; i++){
			

				Move move = chooseMove();
			
					
				/*
				//TESTING			
				if(!ped.sanityCheck()){
					System.out.println(String.format("(%s,%d, %d)", move.name, t, i));
				
					for(int k=0; k< ped.getNActiveNodes(); k++){
						ped.getNode(k).print();
					}
					//System.out.println();
						
					System.out.println(ped.getNActiveNodes());
					ped.printAdjMat();
					System.out.println();
					
					
						
				}
				*/

				

				move.mcmcMove(ped, heat[t]);
				
				iter++;
				counter++;
				
	
			
			}
			
			//record likelihood
			recordLkhd(ped, iter);
			
			//count num success
			//System.out.print(String.format("%d\t%d\t%d\t%f\n", t, moves[0].nAccept, moves[0].nTried, ped.getLogLikelihood()));
			moves[0].nAccept = 0;
			moves[0].nTried = 0;
			
	
			//if the end lkhd did not change much, complete run
			if(counter%100000==0 && Math.abs(initLkhd - ped.getLogLikelihood()) < stopThresh)
				break;
			
			
		}
		
		

		
		
	}
	
	
	
	private void runSample(){
		
		//System.out.println("Sampling...");
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			//record best likelihood
			double currLkhd = ped.getLogLikelihood();
			if(currLkhd > this.bestLkhd){
				this.bestLkhd = currLkhd;
				//System.out.println(bestLkhd);
			}
			
			

			sample(ped);
			//sample2(ped);

			
			//update	
			Move move = chooseMove();
			move.mcmcMove(ped, heat[heat.length-1]);

	
		}
		
		
		//close outfile
		cranefootFamWriter.close();
		convWriter.close();
		//pairWriter.close();
		
	}
	
	

	
	private Move chooseMove(){
		
		double u = rGen.nextDouble();
		double cumProb = 0;
		
		for(Move move : moves){
			
			cumProb += move.getProb();
			
			if(u < cumProb)
				return move;
			
		}
		
		throw new RuntimeException("No move chosen");
		
		
	}
	
	
	
	
	//write relationship to file
	private void sample(Pedigree currPedigree){
		
		//pairwise relationship
		//header for this sample
		
		/*
		//TODO for testing
		pairWriter.write(String.format(">\t%.5f\n", currPedigree.getLogLikelihood()));
		
		for(int i=0; i<currPedigree.numIndiv; i++){
			for(int j=i+1; j<currPedigree.numIndiv; j++){
				
				Path rel = currPedigree.getRelationships()[i][j];
				
				pairWriter.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
			}
		}
		*/
		
		
		//write family relationship
		cranefootFamWriter.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		//currPedigree.clearVisit();
		for(int i=0; i<currPedigree.getNActiveNodes(); i++){
			recordCranefootFam(currPedigree.getNode(i), currPedigree);
		}
		
	
	}
	
	
	private void recordCranefootFam(Node ind, Pedigree currPedigree){
		
		String name = ind.fid + "_" + ind.iid;
		String pa = "0";
		String ma = "0";
		String sampleStatus = ind.sampled ? "000000" : "999999";
		String sex = ind.getSex()==1 ? "1" : "7"; 
		
		//if missing individual and sex not constrained
		//currPedigree.clearVisit();
		//if(currPedigree.sexLocked(ind)==false) sex = "4";
			
		//get parent ids
		for(Node parent : ind.getParents()){
			
			//recordFam(parent);
		
			if(parent.getSex()==0)
				ma = parent.fid + "_" + parent.iid;
			else if(parent.getSex()==1)
				pa = parent.fid + "_" + parent.iid;
			else
				throw new RuntimeException("Parent with unknown sex");
			
		}
		
		//if only one parent is present
		if(ind.getParents().size()==1){
			
			//make missing parent
			int missingParentSex = ind.getParents().get(0).getSex()==1 ? 0 : 1;
			Node missingParent = new Node("missingParent", missingParentCounter+"", missingParentSex, false, -1);
			
			//connect temporarily
			currPedigree.connect(missingParent, ind);
			recordCranefootFam(missingParent, currPedigree);
			currPedigree.disconnect(missingParent, ind);
			
			if(missingParentSex==0) ma = String.format("missingParent_%d", missingParentCounter);
			else pa = String.format("missingParent_%d", missingParentCounter);
			
			missingParentCounter++;
			
		}

		
		
		
		
		//write to file
		cranefootFamWriter.write(String.format("%s\t%s\t%s\t%s\t%s\n", name, pa, ma, sex, sampleStatus));
		
		
	}
	
	
	private void recordLkhd(Pedigree currPedigree, int iter){
		
		convWriter.write(String.format("%d\t%.3f\n", iter, currPedigree.getLogLikelihood()));
		
	}
	
	
	
	public static double acceptanceRatio(double proposalLkhd, double currLkhd, double heat){
		
		//accept if uphill
		if(proposalLkhd > currLkhd){
			return 1d;
		}
		
		//reject if frozen
		if(heat==0) return 0d;
		
		//accept with prob
		double prob = Math.exp(heat * (proposalLkhd - currLkhd));
		
		return prob;
	}
	
	

	
	
}
