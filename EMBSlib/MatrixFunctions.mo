within EMBSlib;
package MatrixFunctions
    extends Modelica.Icons.Package;

    function getTaylorFunction
          input Integer nr;
          input Integer nq;
          input Integer nc;
          input Real[nr,nc] M0;
          input Real[nr,nq,nc] M1;
          input Real[nq] q;
          output Real[nr,nc] M;
  protected
          Real[nr,nc] aux;
    algorithm
          M := zeros(nr,nc);
          for i in 1:nq loop
            aux := M1[:,i,:]*q[i];
            M := M+aux;
          end for;
          M := M0 +M;
          annotation(derivative=EMBSlib.MatrixFunctions.getTaylorFunction_d);
    end getTaylorFunction;

        function getGrMatrix
          input Integer nr;//(for Gr=3 for Ge=nq)
          input Integer nq;//nq
          input Integer nc;//3*nq
          input Real[nr,nc] Min;
          input Real[nq] q;
          output Real[nr, nr] Mout;
  protected
          Real[3,3] aux;
          Integer c=1;
        algorithm
          Mout := zeros(3,3);
          for i in 1:nq loop
            c :=(i - 1)*3;//column index
            aux[:,1] := Min[:,c+1];
            aux[:,2] := Min[:,c+2];
            aux[:,3] := Min[:,c+3];
            Mout := Mout+aux*q[i];
          end for;
        end getGrMatrix;

    function getTaylorFunction_d
          input Integer nr;
          input Integer nq;
          input Integer nc;
          input Real[nr,nc] M0;
          input Real[nr,nq,nc] M1;
          input Real[nq] q;
          input Real[nq] der_q;
          output Real[nr,nc] der_M;
  protected
          Real[nr,nc] aux;
    algorithm
          der_M := zeros(nr,nc);
          for i in 1:nq loop
            aux := M1[:,i,:]*q[i];
            der_M := der_M+aux;
          end for;
          der_M := M0 +der_M;

          annotation(derivative(order=2)=EMBSlib.MatrixFunctions.getTaylorFunction_2d);
    end getTaylorFunction_d;

    function getTaylorFunction_2d
          input Integer nr;
          input Integer nq;
          input Integer nc;
          input Real[nr,nc] M0;
          input Real[nr,nq,nc] M1;
          input Real[nq] q;
          input Real[nq] der_q;
          input Real[nq] der_2_q;
          output Real[nr,nc] der_2_M;
  protected
          Real[nr,nc] aux;
    algorithm
          der_2_M := zeros(nr,nc);
          for i in 1:nq loop
            aux := M1[:,i,:]*q[i];
            der_2_M := der_2_M+aux;
          end for;
          der_2_M := M0 +der_2_M;
    end getTaylorFunction_2d;

        function getGeMatrix
          input Integer nr;//(for Gr=3 for Ge=nq)
          input Integer nq;//nq
          input Integer nc;//3*nq
          input Real[nr,nc] Min;
          input Real[nq] q;
          output Real[nq, 3] Mout;
  protected
          Real[nq,3] aux;
          Integer c=1;
        algorithm
          Mout := zeros(nq,3);
          for i in 1:nq loop
            c :=(i - 1)*3;//column index
            aux[:,1] := Min[:,c+1];
            aux[:,2] := Min[:,c+2];
            aux[:,3] := Min[:,c+3];
            Mout := Mout+aux*q[i];
          end for;
        end getGeMatrix;
end MatrixFunctions;
