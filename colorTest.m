cc=hsv(12);
figure(1); 
hold all;
for i=1:12
    plot([0 1],...
        [0 i],...
        'color',cc(i,:));
end