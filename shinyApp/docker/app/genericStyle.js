[

  {"selector": "node",
     "css": {
        "text-valign": "center",
        "text-halign": "center",
        "border-color": "red",
         "background-color": "#FFFFFF",
        "border-width": 1,
        "label": "data(id)",
        "height": 60,  // defaults
        "width": 60
        }},

   {"selector": "edge",
      "css": {
        "width": "1px",
        "line-color": "blue",
        "target-arrow-shape": "triangle",
        "target-arrow-color": "black",
        "curve-style": "bezier"
        }},

   {"selector": "edge:selected", "css": {
      "line-color": "rgb(200, 200, 0)",
      "overlay-opacity": 0.2,
      "overlay-color": "gray"
      }},

   {"selector":"node:selected", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "border-width": "3px",
       "overlay-opacity": 0.2,
       "overlay-color": "gray"
        }}
   ]
