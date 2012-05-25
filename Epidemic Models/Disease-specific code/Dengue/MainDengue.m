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

% Charger le fichier de données
LoadPath = '/Users/dureaujoseph/Dropbox/Taf/simforence/simforence-population-based/data/dengue/data_forC';
data = load([LoadPath '/data_0.data']);

%Charger le fichier de paramètres
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Dengue';
load([SavePath '/ParametersDengue.mat']) 

%% Définition des structures Data et Dengue_Model
        
Data = struct();
Data.Observations(11,2:length(data)+1) = data; %% l'incidence est calculée dans la 11e composante du vecteur d'état
Parameters.ComputationTStep = 0.05; %% l'unité de temps est le mois. On discrétise dans ce cas le temps à 0.05 mois
Data.Instants = [0 (1:length(data))/Parameters.ComputationTStep];
Data.NbComputingSteps = [0 diff(Data.Instants)];
Data.ObservedVariables = 11*ones(1,length(Data.Instants));
Parameters.DiffusionType = 'Add'; % Ceci indique que la diffusion que suit log(beta_t) est un simple mouvement Brownien

Dengue_Model = struct();
Dengue_Model.EKF_projection = @Dengue_EKF_projection; % fonction qui propage le vecteur d'etat et la covariance entre deux instants d'observation
Dengue_Model.InitializeParameters = @Dengue_Initialize; % fonction qui initialise le vecteur d'etat en fonction de la valeur des paramètres
% Définition des matrices Hk (Observation Jacobian) et Rk
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
%   dont le champs "Estimated" a été fixé à 1 dans le script
%   SaveDengueParameters.m
%   Pendant que le simplexe tourne deux figures sont tracées: les courbes
%   d'incidence observée (vert) et estimée (bleu), aindi que la trajectoire
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




