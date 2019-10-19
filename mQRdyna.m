function [Q, D, ite] = mQRdyna(A, k, tol)
  [Q, Am]=hess(A);
  l=length(Am);
  ite=0; %esto es necesario para resolver ciertos problemas del proy.
  for i=l:-1:2
    Qm=eye(l);
    [Qm(1:i,1:i), Am(1:i,1:i), it]=mp(Am(1:i,1:i), k, tol); %mp = macropaso
    ite = ite + 1;
     %agregar 1 al contador 'ite'.
    Q=Q*Qm; %acumulacion de las Q's.
  end
  lams=diag((Q')*A*Q);
  [v, P] = sort(abs(lams), 'descend'); %ordenar los valores propios.
  D=diag(lams(P));
  Q=Q(:,P); %permutacion.
end

function [Q, Am, i] = mp(A,k, tol)
%Macro-Paso con el shift de Wilkinson, funcionara siempra para matrices
%simetricas.
    Am=A;
    m=1;
    L=length(Am);
    Q=eye(L);
    for i=1:k
      b=Am(L-1,L);
      if abs(b)<tol
        break;
      end
      %Calculos del shift de wilkinson incluidos aqui para hacerlo dinamico
      a=Am(L-1,L);
      c=Am(L,L);
      %como la sub matriz de 2x2 es simetrica, solo necesitamos a y c.
      delta=(a-c)/2;
      rho=c-(sign(delta)*b^2)/(abs(delta)+sqrt(delta^2+b^2)); %wilkinson.
      [Qm,R]=qr(Am-rho*eye(L)); %algoritmo QR.
      Am=R*Qm+rho*eye(L); %Paso iterativo del algoritmo QR, regresando
      %el shift.
      Q = Q*Qm; %acumulacion de Q.
    end
end