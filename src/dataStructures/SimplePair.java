package dataStructures;

public class SimplePair<H, T> {
	private H first;
	private T second;
	
	public SimplePair(H first, T second){
		this.first = first;
		this.second = second;
	}

	public H getFirst() {
		return first;
	}

	public T getSecond() {
		return second;
	}
	
	public void set(H h, T t){
		first = h;
		second = t;
	}
	
	

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((first == null) ? 0 : first.hashCode());
		result = prime * result + ((second == null) ? 0 : second.hashCode());
		return result;
	}

	@SuppressWarnings("unchecked")
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		
		SimplePair<H, T> other = (SimplePair<H, T>) obj;
		if (first == null) {
			if (other.first != null)
				return false;
		} else if (!first.equals(other.first))
			return false;
		if (second == null) {
			if (other.second != null)
				return false;
		} else if (!second.equals(other.second))
			return false;
		return true;
	}
	
	
}
