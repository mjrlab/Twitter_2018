%2017-12-04, EL:
%   - write function to study similarity between two parabolas via
%   cross-correlation.
%   - Fix one curve and move the other along it. Have same number of
%   points in xcorr.

%%
function [delayTime, alignedTestIx, refIx, lagIx, mycorr] = myXcorr(test, ref, interval, maxlag, TOPLOT)
%Compute delayTime between two curves by finding delay that maximizes 
%cross correlation between them. 
%
%Inputs:
%   - test = [x1 y1; x2 y2; ...], ref = [x1 y1; x2 y2; ...]
%   - interval = [lo hi] in units of x's in test; 
%   - maxlag = number (units of x1).
%       - assuming that x's are equally spaced in test and ref
%   - TOPLOT = 1 means plot alignment.
%   
%Output:
%   - delayTime (hrs)
%   - alignedTestIx = [ix1 ix2 ... ixN] indices of points on test curve 
%   that aligns best with the ref curve in [lo hi] interval
%   - refIx = [ix1 ix2 ... ixN] indices of points on ref curve in the [lo
%   hi] interval
%   - lagIx = [ix1 ... ixN] changes to [loIx hiIx] on test curve that have
%   been used to compute correlations; e.g., ix1 = -10 means that
%   test([loIx:hiIx] - 10) was compared against ref(refIx)
%   - mycorr = [corr1 ... corrN] list of correlation coefficients from
%   sliding the correlation window according to lagIx. 

%run in test mode
TOTEST = 0;
if TOTEST == 1
warning('myXcorrFun running in test mode!');   
    
x = -9:0.25:9;
y1 = (x-2).^2;
y2 = (x+2).^2;

test = [x' y2'];
ref = [x' y1'];
interval = [-2 2];
maxlag = 6;

TOPLOT = 1;
end

%ref curve that falls w/in interval
fullref = ref;
ref_lo_ix = find(ref(:,1) > interval(1),1,'first');
ref_hi_ix = find(ref(:,1) < interval(2),1,'last');
refIx = [ref_lo_ix:ref_hi_ix]; %output var
ref = ref(ref_lo_ix:ref_hi_ix,:);
%ref = ref(ref(:,1) > interval(1) & ref(:,1) < interval(2),:);

%endpoints of test curve w/in interval
test_lo_ix = find(test(:,1) > interval(1),1,'first');
test_hi_ix = find(test(:,1) < interval(2),1,'last');

%range w/in which the hi endpoint will vary
max_test_hi_ix = find(test(:,1) < test(test_hi_ix,1) + maxlag, 1, 'last');
min_test_hi_ix = find(test(:,1) > test(test_hi_ix,1) - maxlag, 1, 'first');

%lags in 'ix' space: by how many indices to translate the default 'test'
%window
lag_ix = [min_test_hi_ix:max_test_hi_ix] - test_hi_ix;
lagIx = lag_ix; %output var

%compute correlation between [ref_window] and [test_window]+lag for all
%lags
for l = 1:numel(lag_ix)
    %get endpoints of test interval for the current lag
    lo = test_lo_ix + lag_ix(l);
    hi = test_hi_ix + lag_ix(l);
    
    %grab shifted test region of same length as ref region
    myref = ref;
    mytest = test(lo:hi,:); 
    
    %check lengths
    assert(numel(myref(:,1)) == numel(mytest(:,1)));
    
    %compute correlation coef
    rho = corrcoef(myref(:,2),mytest(:,2)); %returns a 2x2 matrix, grab off-diagonal entry
    
    %rho = sum(myref(:,2).*mytest(:,2))/(std(myref(:,2))*std(mytest(:,2)));
    %out_rho(l) = rho;
    
    %save output
    out_lo(l) = lo;
    out_hi(l) = hi;
    out_rho(l) = rho(1,2);
    
end

%find top lag
[~,m_ix] = max(out_rho);
mycorr = out_rho;

%find endpoints of best aligned interval on test curve (in 'ix' units)
alignedTestIx = lag_ix(m_ix) + [test_lo_ix:test_hi_ix];

%find shift between test curve and best-aligned test curve (in time units)
delayTime = test(test_lo_ix,1) - test(lag_ix(m_ix) + test_lo_ix,1);

%plot to test
if isempty(TOPLOT) 
    TOPLOT = 0;
end
if TOPLOT == 1
figure();
%corr profile
subplot(3,1,1);
plot(lag_ix,out_rho,'bo-');
xlabel('lag');
ylabel('corr');

%original curves, highlight region on ref curve
ax2=subplot(3,1,2);
plot(ref(:,1),ref(:,2),'ks-','DisplayName','ref');
hold on;
plot(fullref(:,1),fullref(:,2),'k-','DisplayName','');
%plot(test(test_lo_ix:test_hi_ix,1),test(test_lo_ix:test_hi_ix,2),'gs-');
plot(test(:,1),test(:,2),'b-','DisplayName','test');
xlabel('time');
ylabel('y');
legend;

%ref curve region used for alignment, region on test curve that best aligns
%with it (also aligned version)
ax3=subplot(3,1,3);
plot(ref(:,1),ref(:,2),'ks-','DisplayName','ref');
hold on;
plot(fullref(:,1),fullref(:,2),'k-','DisplayName','');
plot(test(:,1),test(:,2),'b--','DisplayName','test');
plot(test(:,1)+delayTime,test(:,2),'b.-',...
    'markerfacecolor','b','markersize',4,'DisplayName','aligned test');
plot(test(alignedTestIx,1)+delayTime,test(alignedTestIx,2),'bs',...
    'markerfacecolor','b','markersize',4,'DisplayName','aligned test');
xlabel('time');
ylabel('aligned');
legend;
title(['lag = ' num2str(delayTime,'%2.2f') ' time units']);
linkaxes([ax2,ax3],'xy');
end

end
