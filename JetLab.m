%Keaton Harlow

clear; clc;

%LabData.xlsx needs to placed in the same as this file to run
%Reads in data from the excel spreadsheet
filename = 'LabData.xlsx';
data = xlsread(filename);

%Corrects data so that temperature is in K and pressure is absolute
for i = 1 : 1 : length(data(:, 1))
   for j = 2 : 1 : 6
      data(i, j) = data(i , j) + 273; 
   end
   for j = 7 : 1 : 11
      data(i, j) = data(i, j) + 101.325; 
   end
end

%Sets asumption the the efficiency of the nozzle is .98
eta_nozzle = .98;
%Assuming air as working fluid, therefore R is\
%approximtely 287 kJ kg^-1 K^-1
R = 287;
%Assuming air as working fluid, gamma = 1.4
gamma = 1.4;
%Assuming air as working fluid with no dependence on temperature
%cp = 1004.5 kJ kg^-1 K^-1
cp = 1004.5;
%Conversion factors
LphtoJps = 1 / 3600 / 1000 * 43400000 * 832;
Lphtokgps = 1 / 3600 / 1000 * 832;
%Area of the exit in m^2
A = .0020268299;

%The velocity at the exit: m s^-1
V5 = zeros(length(data(:, 1)), 1);
%Mach at the exit: unitless
M5 = zeros(length(data(:, 1)), 1);
%Static temperature at the exit: Pa
T5 = zeros(length(data(:, 1)), 1);
%Static pressure at the exit: Pa
p5 = zeros(length(data(:, 1)), 1);
%Calculated density, assumption: air as working fluid: kg m^-3
rho5 = zeros(length(data(:, 1)), 1);
%Mass flow rate of air, no fuel: kg s^-1
m_dot = zeros(length(data(:, 1)), 1);
%Internal efficiency: unitless
eta_int = zeros(length(data(:, 1)), 1);
%Theoretical thrust: N
T = zeros(length(data(:, 1)), 1);
%Total mass flow rate including fuel: kg s^-1
m_dot_total = zeros(length(data(:, 1)), 1);
%Heat addition rate: J s^-1
q_dot = zeros(length(data(:, 1)), 1);
%Theoretical stagnation temperature after the combustor
T04 = zeros(length(data(:, 1)), 1);

%Stagnation enthalpy at point 1: kJ kg^-1
h01 = zeros(length(data(:, 1)), 1);
%Stagnation enthalpy at point 2: kJ kg^-1
h02 = zeros(length(data(:, 1)), 1);
%Stagnation enthalpy at point 3: kJ kg^-1
h03 = zeros(length(data(:, 1)), 1);
%Stagnation enthalpy at point 4: kJ kg^-1
h04 = zeros(length(data(:, 1)), 1);
%Isentropic stagnation pressure at point 2: Pa
pr02s = zeros(length(data(:, 1)), 1);
%Isentropic stagnation temperature at point 2: K
T02s = zeros(length(data(:, 1)), 1);
%Isentropic stagnation enthalpy at point 2: kJ kg^-1
h02s = zeros(length(data(:, 1)), 1);
%Compressor efficiency: unitless
eta_c = zeros(length(data(:, 1)), 1);
%Isentropic stagnation pressure at point 4: Pa
pr04s = zeros(length(data(:, 1)), 1);
%Isentropic stagnation temperature at point 4: K
T04s = zeros(length(data(:, 1)), 1);
%Isentropic stagnation enthalpy at point 4: kJ kg^-1
h04s = zeros(length(data(:, 1)), 1);
%Turbine efficiency: unitless
eta_t = zeros(length(data(:, 1)), 1);

for i = 1 : 1 : length(V5)
    %Velocity of the exit
    %Equation: V =
    %sqrt(2*eta_nozzle*gamma/(gamma-1)*R*T4*(1-(p_a/p_4))^((gamma-1)/gamma))
    %real() used due to discrepencies in the data
    V5(i) = real(sqrt(eta_nozzle * data(i, 5) * R * 2 * gamma / (gamma - 1) * (1 - (101.325 / data(i , 10)) ^ ((gamma - 1) / gamma))));
    %Using the T5 and calculated V5 to calculate Mach number at the exit
    M5(i) = V5(i) / sqrt(gamma * R * data(i, 5));
    %Static pressure with gamma = 1.4: p=p0*(1+M5^2/5)^(-7/2)
    p5(i) = data(i, 11) * (1 + M5(i) ^ 2 / 5) ^ (-7 / 2);
    %Static pressure with cp = 1004.5: T=T0-V^2/(2*cp)
    T5(i) = data(i, 6) * (1 + M5(i) ^ 2 / 5) ^ -1;%data(i, 6) - V5(i) ^ 2 / (2 * 1004.5);
    %density assuming ideal gas: rho=p/RT
    rho5(i) = p5(i) * 1000 / (287 * T5(i));
    %mass flow rate: m_dot=rho*A*V
    m_dot(i) = rho5(i) * V5(i) * A;
    %Total mass flow rate: m_dot_total=m_dot_air+m_dot_fuel
    m_dot_total(i) = m_dot(i) + data(i, 12) * Lphtokgps;
    %Heat addition rate
    q_dot(i) = data(i, 12) * LphtoJps;
    %Internal efficiency: eta_int=.5*m_dot*v5^2/(heat addition rate)
    eta_int(i) = .5 * m_dot(i) * V5(i) ^ 2 / q_dot(i);
    %Theoretical thrust: T=m_dot*V
    T(i) = m_dot_total(i) * V5(i);
    %Theoretical temperature at point 4
    %Assumptions: Velocity at point 3 and 4 are equal to 0
    T04(i) =  q_dot(i) / (cp * m_dot(i)) + data(i, 4);
    %Stagnation enthalpy: h0=-69.09*1.1406*T0
    h01(i) = (-69.09*1.1406*data(i,2));
    %Stagnation enthalpy: h0 = -69.09*1.1406*T0
    h02(i) = (-69.09*1.1406*data(i,3));
    %Stagnation enthalpy: h0 = -69.09*1.1406*T0
    h03(i) = (-69.09*1.1406*data(i,4));
    %Stagnation enthalpy: h0 = -69.09*1.1406*T0
    h04(i) = (-69.09*1.1406*data(i,5));
    %Isentropic stagnation pressure: pr02s=p02/p01*pr(t02)
    pr02s(i) = (data(i,8)/data(i,7))*((4.7763*10^-10)*(data(i,2))^(3.8055));
    %Isentropic stagnation temperature: T02s=(pr02s/4.7763E-10)^(1/3.8055)
    T02s(i) = (pr02s(i)/(4.7763*10^-10))^(1/3.8055);
    %Isentropic enthalpy: h02s=-69.09*1.1406*T02s
    h02s(i) = (-69.09*1.1406*T02s(i));
    %Compressor efficiency: eta_c=(h02s-h01)/(h02-h02)
    eta_c(i) = (h02s(i)-h01(i))/(h02(i)-h01(i));
    %Isentropic stagnation pressure: pr04s=p04/p03*pr(T04)
    pr04s(i) = (data(i,10)/data(i,11))*((4.7763*10^-10)*(data(i,4))^(3.8055));
    %Isentropic stagnation temperature: T04s=(pr04s/4.7763E-10)^(1/3.8055)
    T04s(i) = (pr04s(i)/(4.7763*10^-10))^(1/3.8055);
    %Isentropic enthalpy: h04s=-69.09*1.1406*T04s
    h04s(i) = (-69.09+1.1406*T04s(i));
    %Turbine efficiency: eta_t=(h03-h04)/(h03-h04s)
    eta_t(i) = (h03(i)-h04(i))/(h03(i)-h04s(i));
end

%Averages of all temperatures and pressures
avgs = zeros(10, 1);
for i = 1 : 1 : 10
    if i == 5
        for j = 1 : 1 : length(V5)
            avgs(i) = avgs(i) + T5(j);
        end
    elseif i == 10
        for j = 1 : 1 :length(V5)
            avgs(i) = avgs(i) + p5(j);
        end
    else
        for j = 1 : 1 : length(V5)
            avgs(i) = avgs(i) + data(j, i + 1);
        end
    end
%     for j = 1 : 1 : length(V5)
%         avgs(i) = avgs(i) + data(j, i + 1);
%     end
    avgs(i) = avgs(i) / length(V5);
end

%Vector of the specific volumes
v = zeros(6, 1);
press = zeros(6, 1);
temps = zeros(6, 1);
s = zeros(6, 1);
s(1) = 2.54175;
for i = 1 : 1 : 5
    v(i) = R * avgs(i) / avgs(i + 5);
    press(i) = avgs(i + 5);
    temps(i) = avgs(i);
end
for i = 2 : 1 : 5
    s(i) = s(i - 1) + cp * log(temps(i) / temps(i - 1)) - R * log(press(i) / press(i - 1));
end
s(6) = s(1);
v(6) = v(1);
press(6) = press(1);
temps(6) = temps(1);
figure(1)
plot(v, press);
figure(2)
plot(s, temps);