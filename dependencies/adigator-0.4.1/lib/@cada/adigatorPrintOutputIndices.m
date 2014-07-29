function adigatorPrintOutputIndices(x)
%function adigatorPrintOutputIndices(x)
% This module is used to print the derivative mapping to file after all of
% the derivative computations have been printed. The indices and sizes this
% file prints allows the derivative variables which have been printed to
% file to be mapped into their proper Jacobians/Hessians/higher-orders
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR ADIGATORDATA
fid = ADIGATOR.PRINT.FID;
indent = ADIGATOR.PRINT.INDENT;
NDstr = sprintf('%1.0d',ADIGATOR.DERNUMBER);

Dind1 = ['cada',NDstr,'dind1'];
xMrow = x.func.size(1);
xNcol = x.func.size(2);
if isinf(xMrow); xMrow = 1; end
if isinf(xNcol); xNcol = 1; end
[FuncName,DPflag] = cadafuncname(x.id);
if DPflag == 1
  for Vcount = 1:ADIGATOR.NVAROFDIFF
    if ~isempty(x.deriv(Vcount).name)
      VMrow = ADIGATOR.VAROFDIFF(Vcount).size(1);
      VNcol = ADIGATOR.VAROFDIFF(Vcount).size(2);
      if isinf(VMrow); VMrow = 1; end
      if isinf(VNcol); VNcol = 1; end
      DerName = cadadername(FuncName,Vcount,x.id);
      PrevDerName = DerName(1:end-length(ADIGATOR.VAROFDIFF(Vcount).name)-1);
      DerIndex = sub2ind([xMrow*xNcol,VMrow*VNcol],...
        x.deriv(Vcount).nzlocs(:,1),x.deriv(Vcount).nzlocs(:,2));
      nzd = size(x.deriv(Vcount).nzlocs,1);
      if strcmp(PrevDerName(end),'.')
        ADIGATORDATA.INDEXCOUNT = ADIGATORDATA.INDEXCOUNT+1;
        INDEXCOUNT = ADIGATORDATA.INDEXCOUNT;
        IndName = sprintf('Index%1.0d',INDEXCOUNT);
        INDEXNAME = sprintf('Gator%1.0dIndices.Index%1.0d',ADIGATOR.DERNUMBER,INDEXCOUNT);
        if xMrow == 1 && xNcol == 1
          % scalar function
          if VMrow == 1 && VNcol == 1
            % scalar variable
            fprintf(fid,[indent,DerName,'_size = 1;\n']);
            fprintf(fid,[indent,DerName,'_location = 1;\n']);
            Location = 1;
          elseif VMrow == 1 || VNcol == 1
            % vector variable
            fprintf(fid,[indent,DerName,'_size = %1.0d;\n'],VMrow*VNcol);
            Location = zeros(nzd,1);
            Location(:) = DerIndex;
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          else
            % matrix variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d];\n'],VMrow,VNcol);
            Location = zeros(nzd,2);
            [Location(:,1),Location(:,2)] = ind2sub([VMrow,VNcol],DerIndex);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          end
        elseif xMrow == 1 || xNcol == 1
          % vector function
          if VMrow == 1 && VNcol == 1
            % scalar variable
            fprintf(fid,[indent,DerName,'_size = %1.0d;\n'],xMrow*xNcol);
            Location = zeros(nzd,1);
            Location(:) = DerIndex;
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          elseif VMrow == 1 || VNcol == 1
            % vector variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d];\n'],xMrow*xNcol,VMrow*VNcol);
            Location = zeros(nzd,2);
            [Location(:,1),Location(:,2)] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          else
            % matrix variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d,%1.0d];\n'],xMrow*xNcol,VMrow,VNcol);
            Location = zeros(nzd,3);
            [Location(:,1),DerIndex] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
            [Location(:,2),Location(:,3)] = ind2sub([VMrow,VNcol],DerIndex);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          end
        else
          % matrix function
          if VMrow == 1 && VNcol == 1
            % scalar variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d];\n'],xMrow,xNcol);
            Location = zeros(nzd,2);
            [Location(:,1),Location(:,2)] = ind2sub([xMrow,xNcol],DerIndex);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          elseif VMrow == 1 || VNcol == 1
            % vector variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d,%1.0d];\n'],xMrow,xNcol,VMrow*VNcol);
            Location = zeros(nzd,3);
            [DerIndex,Location(:,3)] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
            [Location(:,1),Location(:,2)] = ind2sub([xMrow,xNcol],DerIndex);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,';\n']);
          else
            % matrix variable
            fprintf(fid,[indent,DerName,'_size = [%1.0d,%1.0d,%1.0d,%1.0d];\n'],xMrow,xNcol,VMrow,VNcol);
            Location = zeros(nzd,4);
            [DerIndex,DerIndex2] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
            [Location(:,1),Location(:,2)] = ind2sub([xMrow,xNcol],DerIndex);
            [Location(:,3),Location(:,4)] = ind2sub([VMrow,VNcol],DerIndex2);
            fprintf(fid,[indent,DerName,'_location = ',INDEXNAME,'.',IndName,';\n']);
          end
        end
        ADIGATORDATA.INDICES.(IndName) = Location;
      else
        INDEXNAME = sprintf('Gator%1.0dIndices',ADIGATOR.DERNUMBER);
        % previous deriv vector or scaler (is second or subsequent
        % derivative, so won't have any higher dimension than vector
        if VMrow == 1 && VNcol == 1
          % scalar variable
          fprintf(fid,[indent,DerName,'_size = ',PrevDerName,'_size;\n']);
          Dind1 = cadaindprint(DerIndex);
          fprintf(fid,[indent,DerName,'_location = ',PrevDerName,'_location(',Dind1,',:);\n']);
        elseif VMrow == 1 || VNcol == 1
          % vector variable
          fprintf(fid,[indent,DerName,'_size = [',PrevDerName,'_size,%1.0d];\n'],VMrow*VNcol);
          LocationEnd = zeros(nzd,1);
          [TempDer,LocationEnd(:)] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
          Dind1 = cadaindprint(TempDer);
          ADIGATORDATA.INDEXCOUNT = ADIGATORDATA.INDEXCOUNT+1;
          IndName = sprintf('Index%1.0d',ADIGATORDATA.INDEXCOUNT);
          ADIGATORDATA.INDICES.(IndName) = LocationEnd;
          fprintf(fid,[indent,DerName,'_location = [',PrevDerName,'_location(',Dind1,',:), ',INDEXNAME,'.',IndName,'];\n']);
        else
          % matrix variable
          fprintf(fid,[indent,DerName,'_size = [',PrevDerName,'_size,%1.0d,%1.0d];\n'],VMrow,VNcol);
          LocationEnd = zeros(nzd,2);
          [DerIndex,DerIndex2] = ind2sub([xMrow*xNcol,VMrow*VNcol],DerIndex);
          [LocationEnd(:,1),LocationEnd(:,2)] = ind2sub([VMrow,VNcol],DerIndex2);
          Dind1 = cadaindprint(DerIndex);
          ADIGATORDATA.INDEXCOUNT = ADIGATORDATA.INDEXCOUNT+1;
          IndName = sprintf('Index%1.0d',ADIGATORDATA.INDEXCOUNT);
          ADIGATORDATA.INDICES.(IndName) = LocationEnd;
          fprintf(fid,[indent,DerName,'_location = [',PrevDerName,'_location(',Dind1,',:), ',INDEXNAME,'.',IndName,'];\n']);
        end
      end
    end
  end
end
end
