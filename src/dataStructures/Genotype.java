package dataStructures;

// unordered genotype at a single locus
public class Genotype {

	
	private final String geno;
	
	public Genotype(String geno){
		this.geno = geno;
	}
	
	@Override
	public String toString(){
		return geno;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + (geno.charAt(0) + geno.charAt(1)); //unorder genotype
		return result;
	}
	
	@Override
	public boolean equals(Object obj){
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		
		Genotype other = (Genotype) obj;
		if(other.geno.equals(this.geno))
			return true;
		if(other.geno.charAt(0)==this.geno.charAt(1) && other.geno.charAt(1)==this.geno.charAt(0)) //unordered heterozygotes are the same
			return true;
		
		return false;
	}
	
}
