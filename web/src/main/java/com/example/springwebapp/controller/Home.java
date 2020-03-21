package com.example.springwebapp.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.servlet.ModelAndView;

import javax.servlet.http.HttpSession;

/**
 * Created by Denis on 2/20/2016.
 */

@Controller
public class Home {

    @RequestMapping("/")
    public String home () {
        return "index";
    }

    @RequestMapping("/home-with-session")
    public ModelAndView homeWithSession (HttpSession session) {
        ModelAndView mav = new ModelAndView("index");
        String sid = session.getId();
        mav.addObject("sid", sid);
        return mav;
    }
}