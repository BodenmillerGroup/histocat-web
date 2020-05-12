import RouterComponent from "@/components/RouterComponent.vue";
import Vue from "vue";
import Router from "vue-router";

Vue.use(Router);

export default new Router({
  mode: "history",
  base: process.env.BASE_URL,
  routes: [
    {
      path: "/",
      component: () => import(/* webpackChunkName: "start" */ "@/views/main/Start.vue"),
      children: [
        {
          path: "login",
          component: () => import(/* webpackChunkName: "login" */ "@/views/Login.vue"),
        },
        {
          path: "signup",
          component: () => import(/* webpackChunkName: "login" */ "@/views/SignUp.vue"),
        },
        {
          path: "recover-password",
          component: () => import(/* webpackChunkName: "recover-password" */ "@/views/PasswordRecovery.vue"),
        },
        {
          path: "reset-password",
          component: () => import(/* webpackChunkName: "reset-password" */ "@/views/ResetPassword.vue"),
        },
        {
          path: "profile",
          component: RouterComponent,
          redirect: "profile/view",
          children: [
            {
              path: "view",
              component: () => import(/* webpackChunkName: "profile" */ "@/views/main/profile/UserProfile.vue"),
            },
            {
              path: "edit",
              component: () =>
                import(/* webpackChunkName: "profile-edit" */ "@/views/main/profile/UserProfileEdit.vue"),
            },
            {
              path: "password",
              component: () =>
                import(/* webpackChunkName: "profile-password" */ "@/views/main/profile/UserProfileEditPassword.vue"),
            },
          ],
        },
        {
          path: "main",
          component: () => import(/* webpackChunkName: "main" */ "@/views/main/Main.vue"),
          children: [
            {
              path: "admin",
              component: () => import(/* webpackChunkName: "main-admin" */ "@/views/main/admin/Admin.vue"),
              redirect: "admin/users/all",
              children: [
                {
                  path: "users",
                  redirect: "users/all",
                },
                {
                  path: "users/all",
                  component: () =>
                    import(/* webpackChunkName: "main-admin-users" */ "@/views/main/admin/user/AdminUsers.vue"),
                },
                {
                  path: "users/edit/:id",
                  name: "main-admin-users-edit",
                  component: () =>
                    import(/* webpackChunkName: "main-admin-users-edit" */ "@/views/main/admin/user/EditUser.vue"),
                },
                {
                  path: "users/create",
                  name: "main-admin-users-create",
                  component: () =>
                    import(/* webpackChunkName: "main-admin-users-create" */ "@/views/main/admin/user/CreateUser.vue"),
                },

                {
                  path: "experiments",
                  redirect: "experiments/all",
                },
                {
                  path: "experiments/all",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-admin-experiments" */ "@/views/main/admin/experiment/AdminExperiments.vue"
                    ),
                },
                {
                  path: "experiments/edit/:id",
                  name: "main-admin-experiments-edit",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-admin-experiments-edit" */ "@/views/main/admin/experiment/EditExperiment.vue"
                    ),
                },
              ],
            },
            {
              path: "experiments",
              name: "main-experiments",
              component: () => import(/* webpackChunkName: "main-experiments" */ "@/views/main/ExperimentsView.vue"),
            },
            {
              path: "experiments/create",
              name: "main-experiment-create",
              component: () =>
                import(
                  /* webpackChunkName: "main-experiment-create" */ "@/views/main/admin/experiment/CreateExperiment.vue"
                ),
            },
            {
              path: "experiments/:experimentId/edit",
              name: "main-experiment-edit",
              component: () =>
                import(
                  /* webpackChunkName: "main-experiment-edit" */ "@/views/main/admin/experiment/EditExperiment.vue"
                ),
            },
            {
              path: "experiments/:experimentId/share",
              name: "main-experiment-share",
              component: () =>
                import(/* webpackChunkName: "main-experiment-share" */ "@/views/main/experiment/ShareExperiment.vue"),
            },
            {
              path: "experiments/:experimentId",
              name: "main-experiment",
              component: () =>
                import(/* webpackChunkName: "main-experiment" */ "@/views/main/experiment/ExperimentView.vue"),
              redirect: "experiments/:experimentId/image",
              children: [
                {
                  path: "image",
                  name: "main-experiment-image",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-experiment-image" */ "@/views/main/experiment/image/ImageView.vue"
                    ),
                },
                {
                  path: "data",
                  name: "main-experiment-data",
                  component: () =>
                    import(/* webpackChunkName: "main-experiment-data" */ "@/views/main/experiment/data/DataView.vue"),
                },
              ],
            },
          ],
        },
      ],
    },
    {
      path: "/*",
      redirect: "/",
    },
  ],
});
