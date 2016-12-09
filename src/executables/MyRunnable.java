package executables;

import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import mcmc.SimulatedAnnealing;

public class MyRunnable implements Runnable{
	
	int threadNum;
	Pedigree ped;
	PairwiseLikelihoodCoreStreamPed core;
	SimulatedAnnealing sa;

	
	public MyRunnable(int threadNum, Pedigree ped, PairwiseLikelihoodCoreStreamPed core, SimulatedAnnealing sa){
		
		this.threadNum = threadNum;
		this.ped = ped;
		this.core = core;
		this.sa = sa;

		
	}
	
	
	public void run(){
		
		System.out.println(String.format("Running SA: thread %d", threadNum));
		
		//run MCMC
		//double startTime = System.nanoTime();
		sa.run();
		//double endTime = System.nanoTime();

		//double duration = (endTime - startTime)/1e9; 

		//System.out.println(String.format("Number of singletons: %d", ped.nSingletons[ped.curr]));
		//System.out.println(String.format("Running time: %.1f seconds", duration));
		
		System.out.println(String.format("Finished SA: thread %d", threadNum));

		
	}

}
