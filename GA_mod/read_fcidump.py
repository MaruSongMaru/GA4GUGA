class FCIDUMPReader:
    def __init__(self, filename):
        """Initialize the FCIDUMP reader with a filename."""
        self.filename = filename
        self.two_body_integrals = {}  # (i,j,k,l) -> value
        self.one_body_integrals = {}  # (i,j) -> value
        self.orbital_energies = {}    # i -> value
        self.core_energy = 0.0
        self.read_file()
    
    def read_file(self):
        """Read the FCIDUMP file and process the integrals."""
        with open(self.filename, 'r') as f:
            # Skip the header (everything up to and including &END)
            for line in f:
                line = line.strip()
                if line.endswith("&END"):
                    break
            
            # Process the integral lines
            for line in f:
                parts = line.strip().split()
                if len(parts) != 5:
                    raise ValueError(f"Invalid line format: {line}")
                
                # Parse the line
                value = float(parts[0])
                i, j, k, l = map(int, parts[1:])
                
                # Categorize and store the integral
                self._store_integral(value, i, j, k, l)
    
    def _store_integral(self, value, i, j, k, l):
        """Store integral in the appropriate category based on indices."""
        if i > 0 and j > 0 and k > 0 and l > 0:
            # Two-body integral
            # Store in canonical order to respect 8-fold symmetry
            indices = self._canonical_indices(i, j, k, l)
            self.two_body_integrals[indices] = value
        elif i > 0 and j > 0 and k == 0 and l == 0:
            # One-body integral
            # Store in canonical order (i≤j)
            if i <= j:
                self.one_body_integrals[(i, j)] = value
            else:
                self.one_body_integrals[(j, i)] = value
        elif i > 0 and j == 0 and k == 0 and l == 0:
            # Orbital energy
            self.orbital_energies[i] = value
        elif i == 0 and j == 0 and k == 0 and l == 0:
            # Core energy
            self.core_energy = value
    
    def _canonical_indices(self, i, j, k, l):
        """
        Return canonical ordering of indices to respect 8-fold symmetry:
        ijkl = jilk = ijlk = ... = klij = ...
        
        The canonical ordering is chosen as the one where:
        1. max(i,j) is compared with max(k,l)
        2. The pair with larger maximum comes first
        3. Within each pair, the larger index comes first
        """
        # Sort i,j and k,l pairs internally (larger first)
        ij = (max(i, j), min(i, j))
        kl = (max(k, l), min(k, l))
    
        # Compare the maximum values
        if ij[0] > kl[0]:
            # i,j has larger maximum, so it comes first
            return (ij[0], ij[1], kl[0], kl[1])
        elif kl[0] > ij[0]:
            # k,l has larger maximum, so it comes first
            return (kl[0], kl[1], ij[0], ij[1])
        else:
            # Maxima are equal, compare the minima
            if ij[1] >= kl[1]:
                # i,j has larger or equal minimum, so it comes first
                return (ij[0], ij[1], kl[0], kl[1])
            else:
                # k,l has larger minimum, so it comes first
                return (kl[0], kl[1], ij[0], ij[1])
    
    def get_integral(self, i, j, k, l):
        """
        Get the integral value for given indices.
        For two-body integrals, provide all four indices.
        For one-body integrals, provide i, j and set k=l=0.
        For orbital energies, provide i and set j=k=l=0.
        For core energy, set i=j=k=l=0.
        """
        if i > 0 and j > 0 and k > 0 and l > 0:
            # Two-body integral with 8-fold symmetry
            indices = self._canonical_indices(i, j, k, l)
            return self.two_body_integrals.get(indices, 0.0)
        elif i > 0 and j > 0 and k == 0 and l == 0:
            # One-body integral
            if i <= j:
                return self.one_body_integrals.get((i, j), 0.0)
            else:
                return self.one_body_integrals.get((j, i), 0.0)
        elif i > 0 and j == 0 and k == 0 and l == 0:
            # Orbital energy
            return self.orbital_energies.get(i, 0.0)
        elif i == 0 and j == 0 and k == 0 and l == 0:
            # Core energy
            return self.core_energy
        else:
            raise ValueError("Invalid indices for integral retrieval.")
    
    def summarize(self):
        """Print a summary of the integrals read from the file."""
        print(f"FCIDUMP file: {self.filename}")
        print(f"Core energy: {self.core_energy}")
        print(f"Number of orbital energies: {len(self.orbital_energies)}")
        print(f"Number of one-body integrals: {len(self.one_body_integrals)}")
        print(f"Number of two-body integrals: {len(self.two_body_integrals)}")


# Example usage
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python read_fcidump.py <fcidump_file>")
        sys.exit(1)
    
    fcidump_file = sys.argv[1]
    reader = FCIDUMPReader(fcidump_file)
    reader.summarize()
    
    # Example of accessing integrals
    print("\nSome example integral values:")
    print(f"Core energy: {reader.get_integral(0, 0, 0, 0)}")
    
    if reader.orbital_energies:
        first_orb = next(iter(reader.orbital_energies.keys()))
        print(f"Orbital energy for orbital {first_orb}: {reader.get_integral(first_orb, 0, 0, 0)}")
    
    if reader.one_body_integrals:
        first_one_body = next(iter(reader.one_body_integrals.keys()))
        i, j = first_one_body
        print(f"One-body integral ({i},{j}): {reader.get_integral(i, j, 0, 0)}")
    
    if reader.two_body_integrals:
        first_two_body = next(iter(reader.two_body_integrals.keys()))
        i, j, k, l = first_two_body
        i,j,k,l = 1,2,3,4
        print(f"Two-body integral ({i},{j},{k},{l}): {reader.get_integral(i, j, k, l)}")