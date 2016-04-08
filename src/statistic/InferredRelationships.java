package statistic;

import dataStructures.Path;
import dataStructures.Pedigree;

public class InferredRelationships extends Statistic {

	
	public Object getVal(Pedigree currPedigree) {
		
		Path[][] inferred = currPedigree.getRelationships();
		
		return inferred;
	}
	
	
	
}
