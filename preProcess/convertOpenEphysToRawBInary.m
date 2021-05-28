function ops = convertOpenEphysToRawBInary(ops)

% input .continuous files is in data base temp
% output .bin file in the same folder.

fname       = ops.fbinary;
UTmkdir(fname);

fidout      = fopen(fname, 'w');
if(fidout==-1)
    error(['Could not open file: ',fname])
end



% generates a cell (channel) of structs (blocks) of .continuous files across
% all different directories corresponding to different experimental blocks
ch=load(ops.chanMap);
chans=sort(ch.chanMap);
clear fs
fs=cell(ops.Nchan, 1);

for j = 1:ops.Nchan
    ops.chanMap_KiloRaw(ch.chanMap==chans(j))=j;
    for k = 1:length(ops.root)
        sa(k,1) = dir(fullfile(ops.root{k}, sprintf('*CH%d.continuous', chans(j)) ));
    end
    fs{j} = sa;
end
nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
%
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, 1);

fprintf('Concatenating Open-Ephys data to a single binary file.\n')
tic
ops.nSamplesBlocks=nan(1,nBlocks);
for k = 1:nBlocks
    fprintf('block %d of %d \n', k, nBlocks')
    for j = 1:ops.Nchan
        fid{j} = fopen(fullfile(fs{j}(k).folder, fs{j}(k).name));
        % discard header information
        %fseek(fid{j}, 1024, 0);
        if j==1
           % read header:
           NUM_HEADER_BYTES = 1024;
           fseek(fid{j},0,'bof');
           hdr = fread(fid{j}, NUM_HEADER_BYTES, 'char*1');
           info = getHeader(hdr);
           if isfield(info.header, 'version')
              version = info.header.version;
           else
              version = 0.0;
           end
        else
           % discard header information
           fseek(fid{j}, 1024, 0);
        end

    end
    %
    nsamps = 0;
    nblocks = 0;
    flag = 1;
    while 1
        samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
        for j = 1:ops.Nchan
            
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');
            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');
            nbatches        = floor(numel(rawData)/(nSamples+6));
            
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end
        
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
       
        samples         = samples';
        if(isfield(ops,'common_rejection_mode'))
                switch ops.common_rejection_mode
                    case 'none'
                    case 'mean'
                        samples=samples-repmat(mean(samples),size(samples,1),1);
                    case 'median'
                        samples=samples-repmat(median(samples),size(samples,1),1);
                    otherwise
                        error(['Unknown comon rejection mode ',ops.common_rejection_mode])
                end
        end
            
        %fwrite(fidout, samples, 'int16');
	if nsamps==0
     	    all_samples = samples;
	    block_size = size(samples,2);
	else
	    all_samples = cat(2, all_samples, samples);
        end

        nsamps = nsamps + size(samples,2);
        nblocks = nblocks+1;

        if flag==0
            break;
        end
    end

    %%% Begin LBHB special code
    
    % check if need to remove laser artifacts
    remove_laser_artifact_sec = getparm(ops, 'remove_laser_artifact_sec', 0);
    interp_laser_sec = getparm(ops, 'interp_laser_sec', 0);
    if remove_laser_artifact_sec>0
       % need to get events to figure out when laser was turned on/off,
       % so supply baphy events file
       parmfile = fullfile(ops.runs_root, ops.runs{k});
       LoadMFile(parmfile);
       
       fprintf('removing %.3f sec window around laser on/off\n', ...
          remove_laser_artifact_sec);
       if interp_laser_sec>0
          fprintf('interpolating %.5f sec around laser on/off\n', interp_laser_sec);
       end
      
       spikefs = ops.fs;
       for chan = 1:size(all_samples,1)
          tc = double(all_samples(chan,:))';
          tc_out = remove_opto_artifacts(tc,...
             ops.trial_onsets_{k}, spikefs, exptparams, exptevents, ...
             remove_laser_artifact_sec, interp_laser_sec);
          all_samples(chan,:) = int16(tc_out)';
       end
    end
    
    % check if good_trials are specified
    if isfield(ops,'good_trials') && ~isempty(ops.good_trials{k})
       first_trial=min(ops.good_trials{k});
       last_trial=max(ops.good_trials{k});
       if first_trial==1
          first_bin = 1;
       else
          first_bin=ops.trial_onsets_{k}(first_trial);
       end
       if last_trial<length(ops.trial_onsets_{k})
          last_bin = ops.trial_onsets_{k}(last_trial+1)-1;
       else
          last_bin = nsamps;
       end
       fprintf('keeping good trials %d-%d, reducing samples: %d -> %d\n',...
          first_trial, last_trial, size(all_samples,2), last_bin-first_bin+1);
       
       all_samples = all_samples(:,first_bin:last_bin);
       
    end
    
    %%% End LBHB special code

    for j = 1:nblocks
       e=min((block_size*j),size(all_samples,2));
       samples = all_samples(:,(block_size*(j-1)+1):e);
       fwrite(fidout, samples, 'int16');
    end
    
    ops.nSamplesBlocks(k) = nsamps;
    
    for j = 1:ops.Nchan
       fclose(fid{j}); 
    end
    
end
    
fclose(fidout);
fprintf('Done concatenating, ')
toc
fprintf('\n')

end

function info = getHeader(hdr)
   eval(char(hdr'));
   info.header = header;
end

