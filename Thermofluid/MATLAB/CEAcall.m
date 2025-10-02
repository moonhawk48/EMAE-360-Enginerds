function output = CEAcall(b)
AFR = b(1);
rho = b(2);
temp = b(3);
mN2 = b(4);
mO2 = b(5);
mAr = b(6);
mOctane = b(7);
mHeptane = b(8);
mEthanol = b(9);
P = CEA('problem','uv','equilibrium','f/a',AFR,'eq','rho,kg/m3',rho,'reac', ...
    'ox','N2','wt%',mN2,'t(k)',temp,'ox', ...
    'O2','wt%',mO2,'t(k)',temp, ...
    'ox','Ar','wt%',mAr,'t(k)',temp, ...
    'fu','C8H18,isooctane','wt%',mOctane,'t(k)',temp, ...
    'fu','C7H16,n-heptane','wt%',mHeptane,'t(k)',temp, ...
    'fu','C2H5OH','wt%',mEthanol,'t(k)',temp, ...
    'output','end');
output(1,1) = P.output.temperature;
output(1,2) = P.output.pressure;
output(1,3) = P.output.density;
output(1,4) = P.output.enthalpy;
output(1,5) = P.output.entropy;
output(1,6) = P.output.energy;

end
