function [Gt, rock2D, faults, faces, faultFaces2D] = setUpFaultExample(varargin)
% The following uses fault data (2D corrodinate quantities) to map faults
% to grid faces in the Stø grid model.

    opt.coarse_level = 1;
    opt = merge_options(opt, varargin{:});
    
    fmName = 'Stofm';
    [Gt, rock2D] = getFormationTopGrid(fmName,opt.coarse_level);
    faults = readFaultData();
    [faces, act] = mapFaultsToFaces(Gt, faults, 'refineFactor',2000); 
    fault_active = faults(act);
    faultFaces2D = vertcat(faces{:});


    % Handle parallel faults (which may interfere with simulation due to
    % inadequate grid resolution between two faults when treated as sealing
    % or semi-sealing).
        % possible approaches to handle parallel faults:
        % 1. remove one line, keep other
        % 2. make new, single line inbetween two lines
        % 3. leave both as is, and DONT modify perm inbetween lines
    
    % faces{40} and faces{39} are two parallel fault lines in south, on west
    % and east, respectively.

    % faces{47} and faces{48} are two parallel fault lines in mid-north, on
    % north and south, respectively.

    % Using approach ONE: remove one line, keep other. We remove faces{48},
    % but keep faces{47}. (Removing a line really just means we don't
    % adjust its transmissibility.)
    [~, inx48, ~] = intersect(vertcat(faces{:}), faces{48});
    [~, inx39, ~] = intersect(vertcat(faces{:}), faces{39});
    faultFaces2D([inx48; inx39]) = [];



    % Fault faces on south boundary are "patchy" because grid coordinates
    % are slightly off set from fault coordinates. Yet a fault appears to
    % run along this entire boundary. Thus, we remove fault faces detected
    % near south boundary (i.e., fault line 25) since it does not lie
    % exactly along the boundary and then we add the fault_faces returned
    % by faults_on_boundaried_Sto() to the list of faultFaces2D.
    [~, inx25, ~] = intersect(vertcat(faces{:}), faces{25});
    faultFaces2D([inx25]) = [];
    faultFaces2D = [faultFaces2D; faults_on_boundaries_Sto(Gt)];
    faultFaces2D = unique(faultFaces2D);


end