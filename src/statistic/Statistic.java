package statistic;

import dataStructures.Pedigree;

public abstract class Statistic {

	//take is a Pedigree and return the statistic we're interested in
	public abstract Object getVal(Pedigree currPedigree);

}
