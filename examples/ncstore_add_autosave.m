function ncstore_add_autosave(sim,manager_ref,save_stagename,dt_store)
%NCSTORE_ADD_AUTOSAVE Adds an event to automatically save to an nc file

sim.meta.ncstore_autosave = struct(...
    'manager',manager_ref,...
    'save_stagename',save_stagename,...
    'dt_store',dt_store,...
    'last_store', -dt_store);


sim.append_step_event(@handle_autosave, inf);

end

function handle_autosave(s)
    last_store = s.meta.ncstore_autosave.last_store;
    dt_store = s.meta.ncstore_autosave.dt_store;

    if s.stepping.t > last_store + dt_store - 1e-10
        s.meta.ncstore_autosave.manager.save_time( s, ...
            s.meta.ncstore_autosave.save_stagename,1:length(s.boundaries));
        
        s.meta.ncstore_autosave.last_store = last_store + dt_store * ...
            round((s.stepping.t - last_store)/dt_store);
    end
end