classdef CPRSolverAD < linearSolverAD
    
   methods
        function solver = CPRSolverAD(varargin)
            disp('OK!')
        end

       function result = solveLinearSystem(solver, A, b) %#ok
           error('Not supported')
       end
       
       function [dx, result] = solveLinearProblem(solver, problem)
           
           refSolver = mldivideSolverAD();
           
           % Solve a linearized problem
           isPressure = problem.indexOfPrimaryVariable('pressure');
           
           eqs = problem.equations;
           edd = 1e-2;
           
           % Eliminate the non-cell variables first
           isCell = problem.indexOfType('cell');
           cellIndex = find(isCell);
           cellEqNo = numel(cellIndex);
           
           % Find number of "cell" like variables
           n = problem.getEquationVarNum(cellIndex(1));
           nP = numel(problem);
           
           % Eliminate non-cell variables (well equations etc)
           notCellIndex = find(~isCell);
           
           eliminated = cell(numel(notCellIndex), 1);
           elimNames = problem.equationNames(notCellIndex);
           
           for i = 1:numel(notCellIndex)
               [problem, eliminated{i}] = problem.eliminateVariable(elimNames{i});
           end
           
           
           
           
           diagMod = zeros(n, cellEqNo);
           for i = 1:cellEqNo
               % Find the derivative of current cell block w.r.t the
               % pressure
               pressureJacobi = eqs{cellIndex(i)}.jac{isPressure};
               pressureDiag  = diag(pressureJacobi);
               
               sod = sum(abs(pressureJacobi), 2) - abs(pressureDiag);
                
               diagMod(:,i) = pressureDiag./sod > edd;   
           end
           
           diagMod(all(diagMod == 0, 2), 1) = 1;
           bad = diagMod(:, 1) == 0;
          
           
           % Take first nonzero value instead
           [r,c] = find(diagMod(bad, :));
           firstIndex = accumarray(r,c,[sum(bad), 1], @min);
           diagMod(bad, 1) = firstIndex;
           diagMod = diagMod(:, 1);
           
           % TODO bring up to pairity with cpr
           dx = cell(nP, 1);
           
           
           % Recover non-cell variables
           recovered = false(nP, 1);
           recovered(cellIndex) = true;
           
           % Solve cell equations
           dxCell = refSolver.solveLinearProblem(problem); 
           dx(recovered) = dxCell;
           
           %
           for i = numel(eliminated):-1:1
              pos = notCellIndex(i);
              % OBS: cellEqNo + 1 is a hack WHICH ASSUMES THAT ALL
              % ELIMINATED VARIABLES COMES AFTER CELL VARIABLES
              dVal = recoverVars(eliminated{i}, cellEqNo + 1, dx(recovered));
              dx{pos} = dVal;
              recovered(pos) = true;
           end
           

       end
       
   end
end