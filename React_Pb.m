function R = React_Pb(Pm,Pb,k_adh,Pmax,Rchar)
R = k_adh.*(Pmax - Pb).*Pm./Rchar;