using Mousetrap 
main() do app::Application

    window = Window(app)
    
    # create dummy action
    action = Action("dummy.action", app) do x 
        println("triggered.")
    end

    # model that we will be modifying in the snippet
    root = MenuModel()

    file_submenu = MenuModel()
    add_action!(file_submenu, "Open",action)    
    
    file_recent_submenu = MenuModel()
    add_action!(file_recent_submenu, "Project 01",action)
    add_action!(file_recent_submenu, "Project 02",action)
    add_action!(file_recent_submenu, "Other...",action)
    add_submenu!(file_submenu, "Recent...", file_recent_submenu)
    
    add_action!(file_submenu, "Save Project",action)
    add_action!(file_submenu, "Save Project As...",action)
    
    preparation_submenu = MenuModel()
    calculation_submenu = MenuModel()
    plot_submenu = MenuModel()
    help_submenu = MenuModel()
    
    add_submenu!(root, "File", file_submenu)
    add_submenu!(root, "Preparation", preparation_submenu)
    add_submenu!(root, "Calculation", calculation_submenu)
    add_submenu!(root, "Plots", plot_submenu)
    add_submenu!(root, "Help", help_submenu)
    
    menubar = MenuBar(root)

    # view = PopoverButton(PopoverMenu(root))
    set_margin!(menubar, 10)
    
    center = Label("Welcome to the Potetnial Software!")
    set_margin!(center, 200)
    
    # create a horizontal box
    box = Box(ORIENTATION_HORIZONTAL)
    push_front!(box, center)
    
    # add box to window
    set_child!(window, vbox(menubar, box))

    present!(window)
end