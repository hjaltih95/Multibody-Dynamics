function ycorr = projectPosition(t,ypred)
% project position on constraint surface, EzyRoller
% project the positions on the constraint surface

tol=10^-6;
i=0;
maxiterat = 10;

ycorr= ypred;

[conv,sc,Ds,h] = constraint(ypred);

ap=sc(3:4,:);
while (tol<max(abs(ap)) && i<maxiterat)

    [conv,sc,Ds,h] = constraint(ycorr);
    Dv=Ds(3:4,:);
    ap=sc(3:4,:);
    dy = -Dv.'*((Dv*Dv.')\ap);
    ycorr(1:6) = ycorr(1:6)+dy;
    load=double(100*tol/(abs(ap)));
    i = i+1;
    if i == 10
        error('Max iterations constraints exceded')
    end
end
end
