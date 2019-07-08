<template>
  <div>
    <v-toolbar light>
      <v-toolbar-title>
        Manage Users
      </v-toolbar-title>
      <v-spacer></v-spacer>
      <v-btn color="primary" to="/main/admin/users/create">Create User</v-btn>
    </v-toolbar>
    <v-data-table :headers="headers" :items="users">
      <template slot="items" slot-scope="props">
        <td>{{ props.item.name }}</td>
        <td>{{ props.item.email }}</td>
        <td>{{ props.item.full_name }}</td>
        <td>
          <v-icon v-if="props.item.is_active">mdi-check</v-icon>
        </td>
        <td>
          <v-icon v-if="props.item.is_superuser">mdi-check</v-icon>
        </td>
        <td class="justify-center layout px-0">
          <v-tooltip top>
            <span>Edit</span>
            <v-btn slot="activator" flat :to="{name: 'main-admin-users-edit', params: {id: props.item.id}}">
              <v-icon>mdi-pencil</v-icon>
            </v-btn>
          </v-tooltip>
        </td>
      </template>
    </v-data-table>
  </div>
</template>

<script lang="ts">
  import { userModule } from '@/modules/user';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class AdminUsers extends Vue {
    readonly userContext = userModule.context(this.$store);

    headers = [
      {
        text: 'Name',
        sortable: true,
        value: 'name',
        align: 'left',
      },
      {
        text: 'Email',
        sortable: true,
        value: 'email',
        align: 'left',
      },
      {
        text: 'Full Name',
        sortable: true,
        value: 'full_name',
        align: 'left',
      },
      {
        text: 'Is Active',
        sortable: true,
        value: 'isActive',
        align: 'left',
      },
      {
        text: 'Is Superuser',
        sortable: true,
        value: 'isSuperuser',
        align: 'left',
      },
      {
        text: 'Actions',
        value: 'id',
      },
    ];

    get users() {
      return this.userContext.getters.adminUsers;
    }

    async mounted() {
      await this.userContext.actions.getUsers();
    }
  }
</script>
