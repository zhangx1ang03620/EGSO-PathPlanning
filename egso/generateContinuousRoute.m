
%%
function path=generateContinuousRoute(route,G)

Xmax = size(G,1);
dim = length(route);

path = [route(1) 1];
j=1;
f=1;
failed = false;
while j ~=dim && f<2000 && ~failed
    xa = route(j);
    xb = route(j+1);
    h = abs(xa - xb); % row difference of two adjacent cells
    if h == 0         % same row
        path = [path ; route(j+1) j+1];
    elseif h == 1     % adjacent rows
        path = [path ; route(j+1) j+1];
    else
        k=0;
        while k==0
            if xa <= xb              % cell xa lies above cell xb
                if all(G(xa:xb,j) == 0)    % rows xa..xb of column j are all free
                    path = [path; (xa+1:xb-1)' j*ones(h-1,1);route(j+1) j+1];

                else     % rows xa..xb of column j are not all free
                    w = xa;    % w is the first obstacle cell between xa and xb
                    for k = xa:xb
                        if G(k,j) == 1
                            w = k;
                            break
                        end
                    end
                    if  all(G(w:xb,j+1)==0)    % rows w..xb of column j+1 are all free
                         path = [path; (xa+1:w-1)' j*ones(w-xa-1,1);(w:xb)' (j+1)*ones(xb-w+1,1)];
                    else
                        % ===== Robustness fix: lower bound on the mm search; if no free cell is found,
                        % fail gracefully (path=[]) instead of indexing out of range =====
                        mm=h;
                        found = false;
                        while mm >= 1
                            xb_try = xa+mm-1;
                            if xb_try >= 1 && xb_try <= Xmax && G(xb_try,j+1)==0
                                xb = xb_try;
                                route(j+1)=xb;
                                found = true;
                                break;
                            end
                            mm=mm-1;
                        end
                        if ~found
                            failed = true;
                            break;
                        end
                        j=j-1;
                        continue;
                    end
                end

            else
                if all(G(xa:-1:xb,j) == 0)    % rows xa..xb of column j are all free

                        path = [path; (xa-1:-1:xb+1)' j*ones(h-1,1);route(j+1) j+1];

                else     % rows xa..xb of column j are not all free
                    w = xa;    % w is the first obstacle cell between xa and xb
                    for k = xa:-1:xb
                        if G(k,j) == 1
                            w = k;
                            break
                        end
                    end
                    if  all(G(w:-1:xb,j+1)==0)    % rows w..xb of column j+1 are all free
                         path = [path; (xa-1:-1:w+1)' j*ones(xa-w-1,1);(w:-1:xb)' (j+1)*ones(w-xb+1,1)];
                    else
                        % ===== Robustness fix: same lower bound as above =====
                        mm=h;
                        found = false;
                        while mm >= 1
                            xb_try=xa-mm+1;
                            if xb_try >= 1 && xb_try <= Xmax && G(xb_try,j+1)==0
                                xb = xb_try;
                                route(j+1)=xb;
                                found = true;
                                break;
                            end
                            mm=mm-1;
                        end
                        if ~found
                            failed = true;
                            break;
                        end
                       j=j-1;
                       continue;
                    end
                end
            end
            k=1;
        end
    end
    j=j+1;
    f=f+1;
end
if f>=2000 || failed
    path=[];
end
