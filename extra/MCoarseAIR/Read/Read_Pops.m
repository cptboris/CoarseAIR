%% The Function reads the Level/Group Population Outputted by KONIG
%
%  Input Global Var: - Syst.NAtoms
%                    - Input.FigureFormat ( 'PrePrint' / 'RePrint' / 'Presentation' )
%
function Read_Pops()     

    %%==============================================================================================================
    % 
    % Coarse-Grained method for Quasi-Classical Trajectories (CG-QCT) 
    % 
    % Copyright (C) 2018 Simone Venturi and Bruno Lopez (University of Illinois at Urbana-Champaign). 
    %
    % Based on "VVTC" (Vectorized Variable stepsize Trajectory Code) by David Schwenke (NASA Ames Research Center). 
    % 
    % This program is free software; you can redistribute it and/or modify it under the terms of the 
    % Version 2.1 GNU Lesser General Public License as published by the Free Software Foundation. 
    % 
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    % without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
    % See the GNU Lesser General Public License for more details. 
    % 
    % You should have received a copy of the GNU Lesser General Public License along with this library; 
    % if not, write to the Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA 
    % 
    %---------------------------------------------------------------------------------------------------------------
    %%==============================================================================================================

    global Input Syst Kin Temp
    
    
    for iMol=1:Syst.NMolecules        
        iNBins = Syst.Molecule(iMol).EqNStatesIn;
        
        
        PopFile = strcat(Input.Paths.ToKinRunFldr, '/Pop_', Syst.Molecule(iMol).Name);
        if isfile(strcat(PopFile,'.mat'))
        
            load(strcat(PopFile,'.mat'), 'NSteps', 'PopOverg' )
            Kin.T(Temp.iT).NSteps                  = NSteps;
            Kin.T(Temp.iT).Molecule(iMol).PopOverg = PopOverg; 
            
        else
        

            %% Reading Molecules' Level/Group Populations
            filename = strcat(Input.Paths.ToKinRunFldr, '/pop_', Syst.Molecule(iMol).Name, '.dat')
            delimiter = ' ';
            startRow = 3;
            formatSpec = '%*s%f%*s%*s%[^\n\r]';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            fclose(fileID);
            temp=length(dataArray{:, 1});
            PopVec(1:temp,iMol) = dataArray{:, 1};
            clearvars delimiter startRow formatSpec fileID dataArray ans;


            Kin.T(Temp.iT).Molecule(iMol).PopOverQ = [];
            iStep   = 1;
            iBin    = 1;
            for iVec=1:temp
                if iBin == iNBins+1
                    iStep       = iStep + 1;
                    iBin        = 1;
                else
                    Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,iBin) = PopVec(iVec,iMol);
                    iBin        = iBin + 1;
                end
            end
            Kin.T(Temp.iT).NSteps = iStep-1;
            fprintf('Found %i Steps in the 0-D Solution', Kin.T(Temp.iT).NSteps) 
            
            NSteps   = Kin.T(Temp.iT).NSteps;
            PopOverg = Kin.T(Temp.iT).Molecule(iMol).PopOverg; 
            save(PopFile, 'NSteps', 'PopOverg', '-v7.3');
            
            
        end
        
        
        Kin.T(Temp.iT).Molecule(iMol).Pop = zeros(Kin.T(Temp.iT).NSteps, iNBins);
        for iStep = 1:Kin.T(Temp.iT).NSteps
            Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:) = Kin.T(Temp.iT).Molecule(iMol).PopOverg(iStep,:) .* Syst.Molecule(iMol).T(Temp.iT).GroupsIn.g';
            TotPop   = sum(Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:));
            Kin.T(Temp.iT).Molecule(iMol).DF(iStep,:)  = Kin.T(Temp.iT).Molecule(iMol).Pop(iStep,:) ./ TotPop;
        end
        
        
    end
    
end