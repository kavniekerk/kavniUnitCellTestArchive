Group{
Domain_Dirichlet += Region[750000]; 
Domain_Dirichlet += Region[760000]; 
Domain_Dirichlet += Region[770000]; 
Domain_Dirichlet += Region[780000]; 
Domain_Dirichlet += Region[850000]; 
Domain_Dirichlet += Region[870000]; 
DomainCC_Ele += Region[900000]; 
DomainCC_Ele += Region[910000]; 
DomainCC_Ele += Region[920000]; 
DomainCC_Ele += Region[930000]; 
DomainCC_Ele += Region[960000]; 
DomainCC_Ele += Region[990000]; 
}
Include "MaterialDatabase.pro";
Function {
epsr[Region[900000]] = SteelInd_epsilonr;
epsr[Region[910000]] = SteelInd_epsilonr;
epsr[Region[920000]] = SteelInd_epsilonr;
epsr[Region[930000]] = SteelInd_epsilonr;
epsr[Region[960000]] = Air_epsilonr;
epsr[Region[990000]] = Copper_epsilonr;
}
Constraint { { Name ElectricScalarPotential; Case { { Region Region[750000]; Value -400; } { Region Region[760000]; Value -400; } { Region Region[770000]; Value -400; } { Region Region[780000]; Value -400; } { Region Region[850000]; Value -2000; } { Region Region[870000]; Value 0; } } } }
Constraint { { Name GlobalElectricPotential; Case { } } }
Constraint { { Name GlobalElectricCharge; Case { } } }
Include "Electrostatics.pro";
