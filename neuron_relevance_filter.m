function relevant_neurons = neuron_relevance_filter(peth)
%return boolean array
%relevant_neurons = ones(size(peth,1));
%for num_neuron = 1:size(peth,1)
%    peth_mean = mean(peth(num_neuron,-4000:-2000));
%    peth_stdv = std(peth(num_neuron,-4000:-2000));
%    for mod = (0:step_num)*step_size
%        win_mean = mean(peth(num_neuron,(win_min + mod):(win_min + win_size + mod)));
%        if abs(win_mean - peth_mean) > peth_stdv
%            relevant_neurons(num_neuron) = 0;
%            continue
%        end
%    end
%end