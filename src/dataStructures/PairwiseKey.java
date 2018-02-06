package dataStructures;

//use pair of genotypes as keys; the order of individuals does not matter
//in constructor, the genotypes are ordered so that the smaller allele appears first
public class PairwiseKey {

	
	private String a1;
	private String a2;
	private String b1;
	private String b2;

	
	public PairwiseKey(String a1, String a2, String b1, String b2){
		
		this.a1 = a1;
		this.a2 = a2;
		this.b1 = b1;
		this.b2 = b2;
		
	}
	
	
	public String[] getKey() {
		
		return new String[] {a1, a2, b1, b2};
		
	}
	
	public void setKey(String a1, String a2, String b1, String b2) {
		
		this.a1 = a1;
		this.a2 = a2;
		this.b1 = b1;
		this.b2 = b2;
		
	}
	
	//returns true if this genotype euals other genotype
	private boolean sameGenotypes(String a1_this, String a2_this, String a1_other, String a2_other) {
		
		if(a1_this.equals(a1_other) && a2_this.equals(a2_other))
			return true;
		
		
		else if(a1_this.equals(a2_other) && a2_this.equals(a1_other))
			return true;
		
		return false;
		
		
	}

	
	@Override
	public int hashCode() {
		
		final int prime = 31;
		int result = 1;
		result = prime * result + a1.hashCode();
		result = prime * result + a2.hashCode();
		result = prime * result + b1.hashCode();
		result = prime * result + b2.hashCode();
		return result;
		
	}
	
	@Override
	//keys are equal if unordered genotypes are equal
	public boolean equals(Object obj){
		
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		
		PairwiseKey other = (PairwiseKey) obj;		

		//a's equal each other; b's equal each other
		if(sameGenotypes(this.a1, this.a2, other.a1, other.a2) && sameGenotypes(this.b1, this.b2, other.b1, other.b2))
			return true;
		
		//a equals b; and b equals a
		else if(sameGenotypes(this.a1, this.a2, other.b1, other.b2) && sameGenotypes(this.b1, this.b2, other.a1, other.a2))
			return true;
		
		return false;
		
	}
	
	
	
}
