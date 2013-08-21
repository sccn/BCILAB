% base matrix: access elements: relies on size and mtimes
%
% (c) by Hannes Nickisch, MPI for Biological Cybernetics, 2010 August 23

function Aid = subsref(A,id)

  if ~strcmp(id.type,'()')
    error('Cell contents reference from a non-cell array object.')
  end
  if all(numel(id.subs)~=[1,2])
    error('The object is a 2d matrix.')
  end

  if numel(id.subs)==1            % only a single index => returns a full matrix
    if id.subs{1} == ':'
      Aid = reshape(full(A),[],1);
    else
      m = size(A,1);
      jj = fix(1+(id.subs{1}-1)/m);                 % determine the columns ...
      ii = mod(id.subs{1}-1,m)+1;               % ... and the corresponding rows
      [jju,mm,nn] = unique(jj);                 % those columns will be computed
      Ajju=zeros(size(A,1),numel(jju)); ei=zeros(size(A,2),1);
      for i=1:numel(jju), ei(jju(i))=1; Ajju(:,i)=A*ei; ei(jju(i))=0; end
      Aid = zeros(size(jj));               % allocate memory for return argument
      for i=1:numel(jj), Aid(i) = Ajju(ii(i),nn(i)); end
    end
  else                                  % two indices => returns a matrix object
    ii = 1:size(A,1); if ~all(id.subs{1}==':'), ii = ii(id.subs{1}); end
    jj = 1:size(A,2); if ~all(id.subs{2}==':'), jj = jj(id.subs{2}); end
    Aid = mat([numel(ii),numel(jj)],1,A,ii,jj,'sub');
  end