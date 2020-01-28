package com.example.springwebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;

/**
 * Created by Denis on 2/20/2016.
 */

@Controller
public class Home {

    @RequestMapping("/")
    public String home () {
        return "index";
    }
}