function y = medfilt1nan( x, s )
% standard median filter (copied from medfilt1.m)
n=length(x); 
r=floor(s/2); 
indr=(0:s-1)'; indc=1:n;
ind=indc(ones(1,s),1:n)+indr(:,ones(1,n));
% x0=x(ones(r,1))*0; %Original zero-padding 
x0=NaN(1,r); %New NaN padding
X=[x0'; x'; x0'];
X=reshape(X(ind),s,n); 

%Only fill in if >=3/5 are there.
numnans=sum(isnan(X));
X(:,(numnans>2))=NaN;

%If current is NaN, only fill in if surrounding ones aren't NaN
staynan=isnan(X(3,:)) & (isnan(X(2,:)) | isnan(X(4,:)));
X(:,staynan)=NaN;

%Note other possibilities:
%If either 2 ahead or 2 behind is NaN, just use window of 3 instead of 5??



y=nanmedian(X,1);
end