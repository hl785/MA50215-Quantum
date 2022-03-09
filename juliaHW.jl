# Script like...
print("Hello world!\n")

print("Goodbye world!\n")

# prog like 
mutable struct msgPrinter       # struct = const struct (c++), mutable struct = struct (c++)
    msg::String
    #function msgPrinter(msg)    # constructor, no worries... 
    #    print("Hello ", msg, "!\n")
    #    new(msg)
    #end 
    
    # constructor that defines the destructor? Is destructor seriously a property of the instance instead of the struct!!!
    function msgPrinter(msg)    
       print("Hello ", msg, "!\n")
       x = new(msg)
       f(x) = print("Goodbye ", x.msg, "!\n")   # note f lambda captures x as we never supply x to f when we del x.
       finalizer(f, x)
    end
end

function main()
    x::msgPrinter = msgPrinter("moon")
    x.msg = "sun"
    print(x.msg, "\n")
    finalize(x)
end 

main()