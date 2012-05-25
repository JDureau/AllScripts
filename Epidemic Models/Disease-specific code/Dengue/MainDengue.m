%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Main Dengue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Chargements 

% Charger les fonctions
cd('/Users/dureaujoseph/Dropbox/Taf/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/Dengue'])

% Charger le fichier de donn�es
LoadPath = '/Users/dureaujoseph/Dropbox/Taf/simforence/simforence-population-based/data/dengue/data_forC';
data = load([LoadPath '/data_0.data']);

%Charger le fichier de param�tres
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Dengue';
load([SavePath '/ParametersDengue.mat']) 

%% D�finition des structures Data et Dengue_Model
        
Data = struct();
Data.Observations(11,2:length(data)+1) = data; %% l'incidence est calcul�e dans la 11e composante du vecteur d'�tat
Parameters.ComputationTStep = 0.05; %% l'unit� de temps est le mois. On discr�tise dans ce cas le temps � 0.05 mois
Data.Instants = [0 (1:length(data))/Parameters.ComputationTStep];
Data.NbComputingSteps = [0 diff(Data.Instants)];
Data.ObservedVariables = 11*ones(1,length(Data.Instants));
Parameters.DiffusionType = 'Add'; % Ceci indique que la diffusion que suit log(beta_t) est un simple mouvement Brownien

Dengue_Model = struct();
Dengue_Model.EKF_projection = @Dengue_EKF_projection; % fonction qui propage le vecteur d'etat et la covariance entre deux instants d'observation
Dengue_Model.InitializeParameters = @Dengue_Initialize; % fonction qui initialise le vecteur d'etat en fonction de la valeur des param�tres
% D�finition des matrices Hk (Observation Jacobian) et Rk
% (ObservationMeasurementNoise):
rep1 = Parameters.rep1.Value;
rep2 = Parameters.rep2.Value;
phi = Parameters.phi.Value;
temp = zeros(1,11);
temp(1,11) = rep1*rep2;
Dengue_Model.ObservationJacobian = {};
Dengue_Model.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    Dengue_Model.ObservationJacobian{i} = temp;
end



%%  Simplexe sur Kalman
% 
%   Le simple va optimiser la vraissemblance, en jouant sur les parametres
%   dont le champs "Estimated" a �t� fix� � 1 dans le script
%   SaveDengueParameters.m
%   Pendant que le simplexe tourne deux figures sont trac�es: les courbes
%   d'incidence observ�e (vert) et estim�e (bleu), aindi que la trajectoire
%   de beta_t.
%   
Names = Parameters.Names.Estimated;
ParametersKalman = Parameters;
disp('EKF Maxim. alg. for Hess')
ParametersKalman.Correction = 1;
Initialization = [];
for i = 1:length(Names)
    Initialization(i) = ParametersKalman.(Names{i}).TransfValue ;
end
% ParametersKalman.RWinEKF = 1;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,Dengue_Model,ParametersKalman),Initialization,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',4000));
Names = ParametersKalman.Names.Estimated;
for i = 1:length(Names)
    ParametersKalman.(Names{i}).TransfValue = (x(i));
end
ParametersKalman = UpdateParsTransfToNoTransf(ParametersKalman);
ParametersKalman = Dengue_Model.InitializeParameters(ParametersKalman);
TellParsValues(ParametersKalman)




