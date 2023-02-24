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


%calculate
gamma = 1./sqrt(1-beta.^2);
beta_gamma = gamma.*beta;
dMdbeta = y*M_ej;

deriv_coeff = polyfit(beta, dMdbeta, 2);
final_Mej = M_ej - Mb(deriv_coeff, beta);
figure(1)
plot(beta, dMdbeta, "s", "Color", "#77AC30")
hold on
fit = deriv_coeff(1)*beta.^2 + deriv_coeff(2)*beta + deriv_coeff(3);
Rp2 = calc_R(dMdbeta, fit)
plot(beta, fit, "Color", "#7E2F8E", "LineWidth", 1.5)

legend({"Raw ejecta data", 'Parabolic fit'},'Location','southwest', "FontSize", 7)
title ("$dM_{ej}d\beta$", "Interpreter","latex")
ylabel("$dM_{ej}d\beta [M_\odot]$", "Interpreter","latex")
xlabel("$\beta [c]$", "Interpreter", "latex")

hold off

figure(2)
loglog(beta_gamma, final_Mej, "k*")
hold on

coeff_0 = polyfit(log(beta_gamma), log(final_Mej), 1);
a0 = exp(1)^coeff_0(2)
b0 = coeff_0(1)
M_fit = coeff_0(1)*log(beta_gamma) + coeff_0(2);
loglog(exp(1).^log(beta_gamma), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", "#D95319")
R0 = calc_R(final_Mej, exp(1).^M_fit)

%broken once
[beta_gamma1, beta_gamma2] = splitt2(beta_gamma);
[final_Mej1, final_Mej2] = splitt2(final_Mej);
coeff_broken1 = polyfit(log(beta_gamma1), log(final_Mej1), 1);
a1 = exp(1)^coeff_broken1(2)
b1 = coeff_broken1(1)
M_fit = coeff_broken1(1)*log(beta_gamma1) + coeff_broken1(2);
loglog(exp(1).^log(beta_gamma1), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", 	"#7E2F8E")
R1 = calc_R(final_Mej1, exp(1).^M_fit)

coeff_broken2 = polyfit(log(beta_gamma2), log(final_Mej2), 1);
a2 = exp(1)^coeff_broken2(2)
b2 = coeff_broken2(1)
M_fit = coeff_broken2(1)*log(beta_gamma2) + coeff_broken2(2);
loglog(exp(1).^log(beta_gamma2), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", 	"#7E2F8E")
R2 = calc_R(final_Mej2, exp(1).^M_fit)


%broken twice
[beta_gamma11, beta_gamma12, beta_gamma13] = splitt3(beta_gamma);
[final_Mej11, final_Mej12, final_Mej13] = splitt3(final_Mej);

coeff_broken11 = polyfit(log(beta_gamma11), log(final_Mej11), 1);
a11 = exp(1)^coeff_broken11(2)
b11 = coeff_broken11(1)
M_fit = coeff_broken11(1)*log(beta_gamma11) + coeff_broken11(2);
loglog(exp(1).^log(beta_gamma11), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", "#0072BD")
R11 = calc_R(final_Mej11, exp(1).^M_fit)

coeff_broken12 = polyfit(log(beta_gamma12), log(final_Mej12), 1);
a12 = exp(1)^coeff_broken12(2)
b12 = coeff_broken12(1)
M_fit = coeff_broken12(1)*log(beta_gamma12) + coeff_broken12(2);
loglog(exp(1).^log(beta_gamma12), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", "#0072BD")
R12 = calc_R(final_Mej12, exp(1).^M_fit)

coeff_broken13 = polyfit(log(beta_gamma13), log(final_Mej13), 1);
a13 = exp(1)^coeff_broken13(2)
b13 = coeff_broken13(1)
M_fit = coeff_broken13(1)*log(beta_gamma13) + coeff_broken13(2);
loglog(exp(1).^log(beta_gamma13), exp(1).^M_fit, "--", "LineWidth", 1.5, "Color", "#0072BD")
R13 = calc_R(final_Mej13, exp(1).^M_fit)


title ("$M_{ej}(>\gamma\beta)$", "Interpreter","latex")
ylabel("$M_{ej}(>\gamma\beta) [M_\odot]$", "Interpreter","latex")
xlabel("$\gamma\beta$", "Interpreter", "latex")
legend({"Adapted ejecta data", 'Regular power-law fit', "Broken power-fit with one cutoff", "Broken power-fit with one cutoff", "" + ...
    "Broken power-fit with two cutoffs", "Broken power-fit with two cutoffs", "Broken power-fit with two cutoffs"},'Location','southwest', "FontSize", 7)
set(gca, 'XScale', 'log', 'YScale', 'log');

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
