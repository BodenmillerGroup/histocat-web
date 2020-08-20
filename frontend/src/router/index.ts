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
      component: () => import(/* webpackChunkName: "start" */ "@/views/Start.vue"),
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
          path: "main",
          component: () => import(/* webpackChunkName: "main" */ "@/views/main/Main.vue"),
          children: [
            {
              path: "groups",
              name: "main-groups",
              component: () => import(/* webpackChunkName: "main-groups" */ "@/views/main/group/GroupsView.vue"),
            },
            {
              path: "groups/create",
              name: "main-groups-create",
              component: () =>
                import(/* webpackChunkName: "main-groups-create" */ "@/views/main/group/CreateGroup.vue"),
            },
            {
              path: "groups/:groupId/edit",
              name: "main-groups-edit",
              component: () => import(/* webpackChunkName: "main-groups-edit" */ "@/views/main/group/EditGroup.vue"),
            },
            {
              path: "groups/:groupId",
              name: "main-group",
              component: () => import(/* webpackChunkName: "main-group" */ "@/views/main/group/GroupView.vue"),
              redirect: "groups/:groupId/experiments",
              children: [
                {
                  path: "experiments",
                  name: "main-group-experiments",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-experiments" */ "@/views/main/group/experiment/ExperimentsView.vue"
                    ),
                },
                {
                  path: "experiments/create",
                  name: "main-group-experiments-create",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-experiments-create" */ "@/views/main/group/experiment/CreateExperiment.vue"
                    ),
                },
                {
                  path: "experiments/:experimentId/edit",
                  name: "main-group-experiments-edit",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-experiments-edit" */ "@/views/main/group/experiment/EditExperiment.vue"
                    ),
                },
                {
                  path: "experiments/:experimentId",
                  name: "main-group-experiment",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-experiment" */ "@/views/main/group/experiment/ExperimentView.vue"
                    ),
                },

                {
                  path: "members",
                  name: "main-group-members",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-members" */ "@/views/main/group/members/MembersListView.vue"
                    ),
                },
                {
                  path: "members/create",
                  name: "main-group-members-create",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-members-create" */ "@/views/main/group/members/CreateMember.vue"
                    ),
                },
                {
                  path: "members/:id/edit",
                  name: "main-group-members-edit",
                  component: () =>
                    import(
                      /* webpackChunkName: "main-group-members-edit" */ "@/views/main/group/members/EditMember.vue"
                    ),
                },
              ],
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
                    import(
                      /* webpackChunkName: "profile-password" */ "@/views/main/profile/UserProfileEditPassword.vue"
                    ),
                },
              ],
            },
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
                  path: "users",
                  name: "main-admin-users",
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

                // {
                //   path: "experiments",
                //   redirect: "experiments/all",
                // },
                // {
                //   path: "experiments/all",
                //   component: () =>
                //     import(
                //       /* webpackChunkName: "main-admin-experiments" */ "@/views/main/admin/experiment/AdminExperiments.vue"
                //     ),
                // },
                // {
                //   path: "experiments/edit/:experimentId",
                //   name: "main-admin-experiments-edit",
                //   component: () =>
                //     import(
                //       /* webpackChunkName: "main-admin-experiments-edit" */ "@/views/main/admin/experiment/EditExperiment.vue"
                //     ),
                // },
              ],
            },
            // {
            //   path: "experiments",
            //   name: "main-experiments",
            //   component: () => import(/* webpackChunkName: "main-experiments" */ "@/views/main/ExperimentsView.vue"),
            // },
            // {
            //   path: "experiments/create",
            //   name: "main-experiment-create",
            //   component: () =>
            //     import(
            //       /* webpackChunkName: "main-experiment-create" */ "@/views/main/admin/experiment/CreateExperiment.vue"
            //     ),
            // },
            // {
            //   path: "experiments/:experimentId/edit",
            //   name: "main-experiment-edit",
            //   component: () =>
            //     import(
            //       /* webpackChunkName: "main-experiment-edit" */ "@/views/main/admin/experiment/EditExperiment.vue"
            //     ),
            // },
            // {
            //   path: "experiments/:experimentId/share",
            //   name: "main-experiment-share",
            //   component: () =>
            //     import(/* webpackChunkName: "main-experiment-share" */ "@/views/main/experiment/ShareExperiment.vue"),
            // },
            // {
            //   path: "experiments/:experimentId",
            //   name: "main-experiment",
            //   component: () =>
            //     import(/* webpackChunkName: "main-experiment" */ "@/views/main/experiment/ExperimentView.vue"),
            // },
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
