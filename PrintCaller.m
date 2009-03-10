function PrintCaller
[STACK] = dbstack;
if length(STACK)==1
    fprintf('Use PrintCaller inside a function.\n');
else
    if length(STACK)==2
        STACK(3).file='[command line]';
    end
    fprintf('%s called by %s\n',STACK(2).file,STACK(3).file);
end