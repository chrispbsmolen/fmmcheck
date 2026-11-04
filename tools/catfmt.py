"""Write a scheme in the catalogue's complete file format (symbol matrices A, B, C; the TriadSet; the
closing self-check line), copied from the layout of the served files."""
import bz2
def symmat(name, n, m):
    return '%s:=Matrix(%d, %d, [%s]):' % (name, n, m, ','.join('[' + ','.join('%s_%d_%d' % (name, i + 1, j + 1) for j in range(m)) + ']' for i in range(n)))
def write_full(path, X, Y, Z, triad_strings):
    r = len(triad_strings)
    with bz2.open(path, 'wt') as f:
        f.write(symmat('A', X, Y) + '\n' + symmat('B', Y, Z) + '\n' + symmat('C', Z, X) + '\n')
        f.write('Tensor:=TriadSet([' + ', '.join(triad_strings) + ']):\n')
        f.write('map(expand,A.B-add(LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,1],Tensor)).A)*LinearAlgebra:-Trace(LinearAlgebra:-Transpose(op([1,i,1,2],Tensor)).B)*LinearAlgebra:-Transpose(op([1,i,1,3],Tensor)),i=1..%d));\n' % r)
