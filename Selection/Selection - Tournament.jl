################################################################################################
## Tournament Selection: General Tournament Selection with constant "select worst" probability

abstract TournamentSelection <: SelectionType

type StrictTournSel <: TournamentSelection
	maximizing::Bool
	order::Symbol
	k::Integer
	StrictTournSel(k::Integer; maximizing = true, order = :original) = new(maximizing, order, k)
	StrictTournSel(; maximizing = true, order = :original) 			 = new(maximizing, order, 2)
end

type RelaxedTournSel <: TournamentSelection
	maximizing::Bool
	order::Symbol
	k::Integer
	worst_wins::Vector
	RelaxedTournSel(k::Integer, worst_wins::Real; maximizing = true, order = :original) = new(maximizing, order, k, fill(worst_wins, k - 1))
	RelaxedTournSel(worst_wins::Vector; maximizing = true, order = :original) 					= new(maximizing, order, length(worst_wins) + 1, worst_wins)
	RelaxedTournSel(worst_wins::Real; maximizing = true, order = :original) 						= new(maximizing, order, 2, fill(worst_wins, 1))
end

TournSel(k::Integer, worst_wins::Real, args...) = RelaxedTournSel(k, worst_wins, args...)
TournSel(worst_wins::Real, args...)			 				= RelaxedTournSel(worst_wins, args...)
TournSel(worst_wins::Vector, args...)					 	= RelaxedTournSel(worst_wins, args...)
TournSel(k::Integer, args...)							 			= StrictTournSel(k, args...)
TournSel(args...)										 						= StrictTournSel(args...)

function selectlocations(info::StrictTournSel, fit::Vector, count::Integer; verbose = false)
	better = info.maximizing ? (>) : (<)
	loc = Array(Integer, count)
	alt = rand(1:length(fit), count, info.k)
	if verbose print("fit = \n$fit\n\nalt = \n$alt\n") end
	for i = 1:count
		loc[i] = alt[i,1]
		for j = 2:info.k
			if better(fit[alt[i,j]], fit[loc[i]])
				loc[i] = alt[i,j]
			end
		end
	end

	orderloc!(info, loc)
end

function selectlocations(info::RelaxedTournSel, fit::Vector, count::Integer; verbose = false)
	better = info.maximizing ? (>) : (<)
	worse = info.maximizing ? (<) : (>)
	loc = Array(Integer, count)
	alt = rand(1:length(fit), count, info.k)
	rand_win = rand(count, info.k-1)
	if verbose print("fit = \n$fit\n\nalt = \n$alt\n\nrand_win = \n$rand_win\n") end
	for i = 1:count
		loc[i] = alt[i,1]
		for j = 2:info.k
			if rand_win[i,j-1] < info.worst_wins[j-1]
				if worse(fit[alt[i,j]], fit[loc[i]])
					loc[i] = alt[i,j]
				end
			else
				if better(fit[alt[i,j]], fit[loc[i]])
					loc[i] = alt[i,j]
				end
			end
		end
	end

	orderloc!(info, loc)
end
