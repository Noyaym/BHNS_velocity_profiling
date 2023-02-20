%import data
link = "C:\Noya\המשך מחקר\angular & velocity profiling\APR4-Q3a75.xlsx";
opts = detectImportOptions(link);
opts.SelectedVariableNames = ["beta","y", "Mej"];
T = readtable(link, opts);
beta = T(:,"beta");
beta = table2array(beta);
y = T(:,"y");
y = table2array(y);
M_ej = T(:,"Mej");
M_ej = table2array(M_ej);
M_ej = M_ej(1);

hold on

%calculate
gamma = 1./sqrt(1-beta.^2);
beta_gamma = gamma.*beta;
dMdbeta = y*M_ej;

deriv_coeff = polyfit(beta, dMdbeta, 2);
final_Mej = M_ej - Mb(deriv_coeff, beta);

plot(beta_gamma, final_Mej, "c*")

a = 1;
b = 1;
%not broken
coeff_0 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma, final_Mej)
M_fit = coeff_0(1)*beta_gamma.^coeff_0(2);
plot(beta_gamma, M_fit, "b--")
R = calc_R(final_Mej, M_fit)

%broken once
[beta_gamma1, beta_gamma2] = splitt2(beta_gamma);
[final_Mej1, final_Mej2] = splitt2(final_Mej);
coeff_broken1 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma1, final_Mej1)
M_fit = coeff_broken1(1)*beta_gamma1.^coeff_broken1(2);
plot(beta_gamma1, M_fit, "m--")
R = calc_R(final_Mej1, M_fit)

coeff_broken2 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma2, final_Mej2)
M_fit = coeff_broken2(1)*beta_gamma2.^coeff_broken2(2);
plot(beta_gamma2, M_fit, "y--")
R = calc_R(final_Mej2, M_fit)


%broken twice
[beta_gamma11, beta_gamma12, beta_gamma13] = splitt3(beta_gamma);
[final_Mej11, final_Mej12, final_Mej13] = splitt3(final_Mej);
coeff_broken11 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma11, final_Mej11)
M_fit = coeff_broken11(1)*beta_gamma11.^coeff_broken11(2);
plot(beta_gamma11, M_fit, "g--")
R = calc_R(final_Mej11, M_fit)

coeff_broken12 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma12, final_Mej12)
M_fit = coeff_broken12(1)*beta_gamma12.^coeff_broken12(2);
plot(beta_gamma12, M_fit, "r--")
R = calc_R(final_Mej12, M_fit)

coeff_broken13 = lsqcurvefit(@(ab,bg)pl(ab,bg),[a, b],beta_gamma13, final_Mej13)
M_fit = coeff_broken13(1)*beta_gamma13.^coeff_broken13(2);
plot(beta_gamma13, M_fit, "k--")
R = calc_R(final_Mej13, M_fit)

legend({"Ejecta velocity measurements", 'not broken', "broken once p1", "broken once p2", "broken twice p1", "broken twice p2", "broken twice p3"},'Location','northeast', "FontSize", 7)

hold off


%functions
function M = Mb(deriv_coeff, beta) 
M = deriv_coeff(1)*beta.^3./3 + deriv_coeff(2)*beta.^2./2 + deriv_coeff(3)*beta;
end

function pl = pl(ab, bg)
pl = ab(1)*bg.^ab(2);
end

function [s1, s2] = splitt2(array)
lx = (length(array));
half = ceil(lx/2);
s1 = array(1:half);
s2 = array(half + 1 : end);
end

function [s1, s2, s3] = splitt3(array)
lx = (length(array));
third = ceil(lx/3);
s1 = array(1:lx-2*third);
s2 = array(lx-2*third+1 : lx-2*third+1+third);
s3 = array(lx-2*third+2+third : end);
end

function Rsq = calc_R(ydata, fundata)
Rsq = 1 - E(ydata, fundata)/vari(ydata);
end

function e = E(ydata, fundata)
e = sum((fundata-ydata).^2);
end

function vrnc = vari(data)
vrnc = sum((data-mean(data)).^2);
end





