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
          component: () => import(/* webpackChunkName: "signup" */ "@/views/SignUp.vue"),
        },
        {
          path: "password-recovery",
          component: () => import(/* webpackChunkName: "password-recovery" */ "@/views/PasswordRecovery.vue"),
        },
        {
          path: "reset-password",
          component: () => import(/* webpackChunkName: "reset-password" */ "@/views/ResetPassword.vue"),
        },

        {
          path: "main",
          component: () => import(/* webpackChunkName: "main" */ "@/views/Main.vue"),
          children: [
            {
              path: "groups",
              name: "groups",
              component: () => import(/* webpackChunkName: "groups" */ "@/views/group/GroupsView.vue"),
            },
            {
              path: "groups/create",
              name: "groups-create",
              component: () => import(/* webpackChunkName: "groups-create" */ "@/views/group/CreateGroup.vue"),
            },
            {
              path: "groups/:groupId/edit",
              name: "groups-edit",
              component: () => import(/* webpackChunkName: "groups-edit" */ "@/views/group/EditGroup.vue"),
            },
            {
              path: "groups/:groupId",
              name: "group",
              component: () => import(/* webpackChunkName: "group" */ "@/views/group/GroupView.vue"),
              redirect: "groups/:groupId/experiments",
              children: [
                {
                  path: "experiments",
                  name: "group-experiments",
                  component: () =>
                    import(/* webpackChunkName: "group-experiments" */ "@/views/group/experiment/ExperimentsView.vue"),
                },
                {
                  path: "experiments/create",
                  name: "group-experiments-create",
                  component: () =>
                    import(
                      /* webpackChunkName: "group-experiments-create" */ "@/views/group/experiment/CreateExperiment.vue"
                    ),
                },
                {
                  path: "experiments/:experimentId/edit",
                  name: "group-experiments-edit",
                  component: () =>
                    import(
                      /* webpackChunkName: "group-experiments-edit" */ "@/views/group/experiment/EditExperiment.vue"
                    ),
                },
                {
                  path: "experiments/:experimentId",
                  name: "group-experiment",
                  component: () =>
                    import(/* webpackChunkName: "group-experiment" */ "@/views/group/experiment/ExperimentView.vue"),
                },

                {
                  path: "members",
                  name: "group-members",
                  component: () =>
                    import(/* webpackChunkName: "group-members" */ "@/views/group/members/MembersListView.vue"),
                },
                {
                  path: "members/create",
                  name: "group-members-create",
                  component: () =>
                    import(/* webpackChunkName: "group-members-create" */ "@/views/group/members/CreateMember.vue"),
                },
                {
                  path: "members/:id/edit",
                  name: "group-members-edit",
                  component: () =>
                    import(/* webpackChunkName: "group-members-edit" */ "@/views/group/members/EditMember.vue"),
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
                  component: () => import(/* webpackChunkName: "profile" */ "@/views/profile/UserProfile.vue"),
                },
                {
                  path: "edit",
                  component: () => import(/* webpackChunkName: "profile-edit" */ "@/views/profile/UserProfileEdit.vue"),
                },
                {
                  path: "password",
                  component: () =>
                    import(/* webpackChunkName: "profile-password" */ "@/views/profile/UserProfileEditPassword.vue"),
                },
              ],
            },

            {
              path: "admin",
              component: () => import(/* webpackChunkName: "admin" */ "@/views/admin/Admin.vue"),
              redirect: "admin/users",
              children: [
                {
                  path: "users",
                  redirect: "users",
                },
                {
                  path: "users",
                  name: "admin-users",
                  component: () => import(/* webpackChunkName: "admin-users" */ "@/views/admin/user/AdminUsers.vue"),
                },
                {
                  path: "users/edit/:id",
                  name: "admin-users-edit",
                  component: () => import(/* webpackChunkName: "admin-users-edit" */ "@/views/admin/user/EditUser.vue"),
                },
                {
                  path: "users/create",
                  name: "admin-users-create",
                  component: () =>
                    import(/* webpackChunkName: "admin-users-create" */ "@/views/admin/user/CreateUser.vue"),
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
